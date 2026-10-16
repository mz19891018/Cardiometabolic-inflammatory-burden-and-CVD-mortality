.libPaths(c('work/Rlib',.libPaths())); library(survey); library(survival); library(foreign)
obj<-readRDS('work/audit_models.rds'); d<-obj$data
options(survey.lonely.psu='adjust')
base<-'age10 + factor(riagendr) + indfmpir + smoking_cat + education + race + pa_correct + bp_correct + dm_correct + lipid_med'
ok<-d$weight>0&d$eligstat==1&is.finite(d$cmibi_sd)&d$time>0
ds<-svydesign(ids=~psu,strata=~stratum,weights=~weight,data=d[d$weight>0&d$eligstat==1,],nest=TRUE)
# Comparable cohort for component-removal sensitivity, fixed full score sample and scale.
q<-quantile(d$cmibi[ok],c(0,.25,.5,.75,1)); d$quartile<-cut(d$cmibi,q,include.lowest=TRUE,labels=c('Q1','Q2','Q3','Q4'))
d$quartile_score<-as.numeric(d$quartile)
ds<-update(ds,quartile=d$quartile[d$weight>0&d$eligstat==1],quartile_score=d$quartile_score[d$weight>0&d$eligstat==1])
f<-as.formula(paste('Surv(time,cvd_death)~quartile+',base))
qm<-svycoxph(f,design=subset(ds,seqn%in%d$seqn[ok]),ties='breslow',model=TRUE,x=TRUE,y=TRUE)
capture.output(summary(qm),file='outputs/quartile_results.txt')
# Survey estimate of component-complete sample size and descriptive means.
dd<-subset(ds,seqn%in%d$seqn[ok]); write.csv(as.data.frame(svymean(~ridageyr+cmibi+crp_mgl+egfr2021,dd,na.rm=TRUE)),'outputs/descriptive_means.csv')
# Inverse probability of complete covariates: explicitly a sensitivity, not multiple imputation.
vars<-all.vars(as.formula(paste('~',base))); d$response<-as.numeric(complete.cases(d[,vars]))
rr<-d[ok,]; rr$income_missing<-as.numeric(is.na(rr$indfmpir)); rr$smoking_missing<-as.numeric(is.na(rr$smoking))
response_model<-glm(response~splines::ns(ridageyr,3)+factor(riagendr)+factor(ridreth1)+factor(cycle)+cmibi_sd+I(cmibi_sd^2),data=rr,weights=weight/mean(weight),family=quasibinomial())
rr$prob_complete<-predict(response_model,type='response'); rr$weight_ipw<-rr$weight/pmax(rr$prob_complete,.05)
cc<-rr[rr$response==1,]; descc<-svydesign(ids=~psu,strata=~stratum,weights=~weight_ipw,data=cc,nest=TRUE)
f<-as.formula(paste('Surv(time,cvd_death)~cmibi_sd+',base))
mm<-svycoxph(f,design=descc,ties='breslow',model=TRUE,x=TRUE,y=TRUE)
b<-coef(mm)['cmibi_sd']; se<-sqrt(vcov(mm)['cmibi_sd','cmibi_sd']); df<-degf(descc)
write.csv(data.frame(n=nrow(cc),events=sum(cc$cvd_death),HR=exp(b),lower=exp(b-qt(.975,df)*se),upper=exp(b+qt(.975,df)*se),p=2*pt(-abs(b/se),df),min_response=min(rr$prob_complete),max_response=max(rr$prob_complete)),'outputs/missingness_ipw_sensitivity.csv',row.names=FALSE)
# Raw fasting adult denominator before the legacy inner merges.
dirs<-c(D='2005-2006年',E='2007-2008年',F='2009-2010年',I='2015-2016年',J='2017-2018年')
den<-lapply(names(dirs),function(s){demo<-read.xport(file.path('E:/NHANES',dirs[s],paste0('DEMO_',s,'.xpt'))); glu<-read.xport(file.path('E:/NHANES',dirs[s],paste0('GLU_',s,'.xpt'))); x<-merge(demo,glu,by='SEQN'); x<-x[x$RIDAGEYR>=20&!is.na(x$WTSAF2YR)&x$WTSAF2YR>0,]; pregnant<-!is.na(x$RIDEXPRG)&x$RIDEXPRG==1; data.frame(cycle=s,raw_fasting_adults=nrow(x),pregnant=sum(pregnant),nonpregnant=sum(!pregnant),legacy_rows=sum(d$cycle==s),missing_after_legacy_inner_merges=sum(!pregnant & !x$SEQN%in%d$seqn))})
write.csv(do.call(rbind,den),'outputs/raw_denominator_check.csv',row.names=FALSE)
r<-read.csv('outputs/cox_results.csv'); take<-match(c('cvd_death_M3','cvd_death_M4','CVD_exclude_early_2years','CVD_no_baseline_CVD_cancer','CVD_without_CRP','CVD_without_eGFR'),r$analysis); r<-r[take,]
png('outputs/cvd_sensitivity_forest.png',width=1800,height=1050,res=180)
par(mar=c(5,15,3,2)); yy<-rev(seq_len(nrow(r)))
plot(r$HR,yy,xlim=c(.7,2.1),ylim=c(.5,nrow(r)+.5),log='x',pch=19,yaxt='n',xlab='Hazard ratio per SD (95% CI)',ylab='',main='CMIBI and cardiovascular mortality')
segments(r$lower,yy,r$upper,yy,lwd=2); abline(v=1,lty=2,col='gray50'); axis(2,at=yy,labels=c('Corrected primary model','Add survey cycle','Exclude deaths within 2 years','Exclude baseline CVD/cancer','Remove CRP','Remove eGFR'),las=1); dev.off()
