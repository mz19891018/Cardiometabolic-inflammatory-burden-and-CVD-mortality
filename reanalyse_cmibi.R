args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) args[1] else 'E:/NHANES'
out <- if(length(args)>1) args[2] else 'outputs'
.libPaths(c('work/Rlib', .libPaths()))
library(survey); library(survival); library(foreign)
options(survey.lonely.psu='adjust')
dir.create(out,showWarnings=FALSE,recursive=TRUE)
sink(file.path(out,'analysis_log.txt'),split=TRUE)
cat('Reanalysis ',as.character(Sys.time()),'\n'); print(sessionInfo())
readx <- function(s, file) {
  dirs <- c(D='2005-2006年',E='2007-2008年',F='2009-2010年',I='2015-2016年',J='2017-2018年')
  x <- read.xport(file.path(root,dirs[s],paste0(file,'_',s,'.xpt')))
  names(x)<-tolower(names(x)); stopifnot(!anyDuplicated(x$seqn)); x
}
d <- read.csv(file.path(root,'analysis/analysis_data.csv'),check.names=FALSE)
stopifnot(!anyDuplicated(d$seqn))
cat('Existing derived adult nonpregnant rows:',nrow(d),'\n')
getv <- function(x,n) if(n %in% names(x)) x[[n]] else rep(NA_real_,nrow(x))
yn <- function(x) ifelse(x==1,1,ifelse(x==2,0,NA_real_))
or01 <- function(a,b) ifelse((!is.na(a)&a==1)|(!is.na(b)&b==1),1,ifelse(!is.na(a)&!is.na(b),0,NA_real_))
rawcheck <- list()
for(nm in c('bp_correct','dm_correct','pa_correct','crp_mgl')) d[[nm]]<-NA_real_
for(s in c('D','E','F','I','J')) {
  ix<-which(d$cycle==s); bp<-readx(s,'BPQ'); di<-readx(s,'DIQ'); pa<-readx(s,'PAQ')
  b<-bp[match(d$seqn[ix],bp$seqn),]; q<-di[match(d$seqn[ix],di$seqn),]; p<-pa[match(d$seqn[ix],pa$seqn),]
  cat('DIQ',s,':',paste(grep('di[q|d]0[578]',names(di),value=TRUE),collapse=','),'\n')
  bm<-yn(getv(b,'bpq050a')); bm[is.na(bm)&(getv(b,'bpq020')==2|getv(b,'bpq040a')==2)]<-0
  insulin<-yn(getv(q,'diq050'))
  pillname<-if('did070'%in%names(q)) 'did070' else 'diq070'
  pills<-yn(getv(q,pillname))
  # Structural skips only: explicit negative diagnosis and prediabetes history.
  skip<-getv(q,'diq010')==2 & (is.na(getv(q,'diq160'))|getv(q,'diq160')==2)
  pills[is.na(pills)&!is.na(skip)&skip]<-0
  dm<-or01(insulin,pills)
  a<-yn(getv(p,if(s=='D') 'pad200' else 'paq650'))
  z<-yn(getv(p,if(s=='D') 'pad320' else 'paq665'))
  d$bp_correct[ix]<-bm; d$dm_correct[ix]<-dm; d$pa_correct[ix]<-or01(a,z)
  cr<-readx(s,if(s%in%c('I','J')) 'HSCRP' else 'CRP')
  crv<-getv(cr,if(s%in%c('I','J')) 'lbxhscrp' else 'lbxcrp')
  if(!s%in%c('I','J')) crv<-crv*10
  d$crp_mgl[ix]<-crv[match(d$seqn[ix],cr$seqn)]
  rawcheck[[s]]<-data.frame(cycle=s,n=length(ix),bp_changed=sum(d$bp_med[ix]!=bm,na.rm=TRUE),dm_changed=sum(d$dm_med[ix]!=dm,na.rm=TRUE),pa_changed=sum(d$phys_activity[ix]!=d$pa_correct[ix],na.rm=TRUE),crp_max_difference=max(abs(d$lbxcrp[ix]*10-d$crp_mgl[ix]),na.rm=TRUE))
}
write.csv(do.call(rbind,rawcheck),file.path(out,'raw_variable_checks.csv'),row.names=FALSE)
d$sbp[d$sbp<=50]<-NA; d$dbp[d$dbp<20]<-NA
components <- function(x,egfr=x$egfr) data.frame(HbA1c=x$lbxgh,TG=log(pmax(x$lbxtr,1)),HDL=-x$lbdhdd,TC=x$lbxtc,SBP=x$sbp,DBP=x$dbp,BMI=x$bmxbmi,waist=x$bmxwaist,CRP=log(pmax(x$lbxcrp,.01)),eGFR=-egfr,urate=x$lbxsua,albumin=-x$lbxsal,ALP=log(pmax(x$lbxsapsi,1)))
standard <- function(v) (v-mean(v,na.rm=TRUE))/sd(v,na.rm=TRUE)
oldcomp<-as.data.frame(lapply(components(d),standard)); d$cmibi_old<-rowSums(oldcomp)
oldok<-complete.cases(oldcomp); d$cmibi_old_sd<-(d$cmibi_old-mean(d$cmibi_old[oldok]))/sd(d$cmibi_old[oldok])
d$time_old<-d$followup_years; d$time<-d$permth_exm/12
d$weight<-d$wt_fast
d$age10<-d$ridageyr/10; d$race<-factor(d$ridreth1)
d$smoke_old<-d$smoking; d$educ_old<-ifelse(d$dmdeduc2%in%c(7,9),NA,d$dmdeduc2)
d$alcohol_old<-ifelse(d$alcohol==9,NA,d$alcohol)
d$psu<-interaction(d$cycle,d$sdmvstra,d$sdmvpsu,drop=TRUE)
d$stratum<-interaction(d$cycle,d$sdmvstra,drop=TRUE)
des <- svydesign(ids=~psu,strata=~stratum,weights=~weight,data=d[is.finite(d$weight)&d$weight>0,],nest=TRUE)
results<-list(); models<-list()
fit <- function(label,event,covs,old=FALSE,select=rep(TRUE,nrow(d))) {
  tm<-if(old) 'time_old' else 'time'
  vars<-unique(c(tm,event,all.vars(as.formula(paste('~',covs)))))
  valid<-complete.cases(d[,vars]) & d[[tm]]>0 & select & d$eligstat==1
  valid[is.na(valid)]<-FALSE
  ids<-d$seqn[valid]; ds<-subset(des,seqn%in%ids)
  f<-as.formula(paste0('Surv(',tm,',',event,') ~ ',covs))
  m<-svycoxph(f,design=ds,ties='breslow',model=TRUE,x=TRUE,y=TRUE)
  nm<-grep('cmibi|score_',names(coef(m)),value=TRUE)[1]
  b<-coef(m)[nm]; se<-sqrt(vcov(m)[nm,nm]); dfree<-degf(ds); crit<-qt(.975,dfree)
  rr<-data.frame(analysis=label,event=event,n=nrow(ds$variables),events=sum(ds$variables[[event]]),df=dfree,HR=exp(b),lower=exp(b-crit*se),upper=exp(b+crit*se),p=2*pt(-abs(b/se),dfree),beta=b,se=se)
  results[[length(results)+1]]<<-rr; models[[label]]<<-m; print(rr); m
}
oldcov<-'cmibi_old_sd + age10 + riagendr + indfmpir + smoke_old + alcohol_old + phys_activity + educ_old + race + bp_med + dm_med + lipid_med'
fit('legacy_covariates_standard_survey','cvd_death',oldcov,TRUE)
# Correct creatinine calibration and use published race-free CKD-EPI 2021.
d$scr_cal<-d$lbxscr; d$scr_cal[d$cycle=='D']<--.016+.978*d$scr_cal[d$cycle=='D']
j<-d$cycle=='J'
d$scr_cal[j]<-1.051*d$lbxscr[j]-.06945
k<-ifelse(d$riagendr==1,.9,.7); a<-ifelse(d$riagendr==1,-.302,-.241)
d$egfr2021<-142*pmin(d$scr_cal/k,1)^a*pmax(d$scr_cal/k,1)^(-1.2)*.9938^d$ridageyr*ifelse(d$riagendr==2,1.012,1)
newraw<-components(d,d$egfr2021); newraw$CRP<-log(pmax(d$crp_mgl,.1))
newraw$albumin[j]<--(1.044*d$lbxsal[j]+.01128)
newraw$ALP[j]<-log(10^(1.001*log10(d$lbxsapsi[j])-.04294))
newraw$urate[j]<-.9323*d$lbxsua[j]+.2326
ref<-is.finite(d$weight)&d$weight>0&d$eligstat==1
params<-data.frame(component=names(newraw),mean=sapply(newraw,function(x)mean(x[ref],na.rm=TRUE)),sd=sapply(newraw,function(x)sd(x[ref],na.rm=TRUE)))
write.csv(params,file.path(out,'reference_parameters.csv'),row.names=FALSE)
nc<-as.data.frame(Map(function(x,mu,sig)(x-mu)/sig,newraw,params$mean,params$sd))
d$cmibi<-rowSums(nc); ok<-ref&complete.cases(nc)
cmu<-mean(d$cmibi[ok]); csd<-sd(d$cmibi[ok]); d$cmibi_sd<-(d$cmibi-cmu)/csd
write.csv(data.frame(mean=cmu,sd=csd),file.path(out,'score_scale.csv'),row.names=FALSE)
d$smoking_cat<-factor(d$smoking,levels=0:2)
d$education<-factor(ifelse(d$dmdeduc2%in%1:2,1,ifelse(d$dmdeduc2==3,2,ifelse(d$dmdeduc2%in%4:5,3,NA))))
d$diabetes_correct<-or01(ifelse(d$diq010==1,1,ifelse(d$diq010%in%c(2,3),0,NA)),d$dm_correct)
d$cycle<-factor(d$cycle)
base<-'age10 + factor(riagendr) + indfmpir + smoking_cat + education + race + pa_correct'
m3<-paste(base,'+ bp_correct + dm_correct + lipid_med')
# Alcohol omitted in primary corrected analysis: ALQ101 and ALQ111 have different reference periods.
des<-svydesign(ids=~psu,strata=~stratum,weights=~weight,data=d[ref,],nest=TRUE)
for(ev in c('mortstat','cvd_death','cancer_death')) {
  fit(paste0(ev,'_M1'),ev,'cmibi_sd + age10 + factor(riagendr)')
  fit(paste0(ev,'_M2'),ev,paste('cmibi_sd +',base))
  fit(paste0(ev,'_M3'),ev,paste('cmibi_sd +',m3))
  fit(paste0(ev,'_M4'),ev,paste('cmibi_sd +',m3,'+ cycle'))
}
fit('CVD_age_nonlinear','cvd_death',paste('cmibi_sd + splines::ns(ridageyr,df=4) +',sub('age10 + ','',m3,fixed=TRUE)))
fit('CVD_legacy_alcohol_sensitivity','cvd_death',paste('cmibi_sd +',m3,'+ factor(alcohol_old)'))
fit('CVD_exclude_early_2years','cvd_death',paste('cmibi_sd +',m3),select=!(d$mortstat==1&d$time<=2))
cv<-as.matrix(d[,paste0('mcq160',letters[2:6])]); d$cvd_strict<-ifelse(rowSums(cv==1,na.rm=TRUE)>0,1,ifelse(rowSums(cv==2,na.rm=TRUE)==5,0,NA))
fit('CVD_no_baseline_CVD_cancer','cvd_death',paste('cmibi_sd +',m3),select=d$cvd_strict==0&d$mcq220==2)
for(drop in c('CRP','eGFR')) {
  v<-rowSums(nc[,setdiff(names(nc),drop)]); d$score_drop<-standard(v)
  des<-update(des,score_drop=d$score_drop[ref])
  fit(paste0('CVD_without_',drop),'cvd_death',paste('score_drop +',m3))
}
# Restricted cubic spline equivalent: natural cubic spline with four fixed quantile knots, 3 df.
kn<-quantile(d$cmibi_sd[ok],c(.05,.35,.65,.95))
bs<-splines::ns(d$cmibi_sd,knots=kn[2:3],Boundary.knots=kn[c(1,4)])
# Orthogonalize nonlinear columns against intercept and linear score, preserving linear column.
ii<-is.finite(d$cmibi_sd); residbasis<-matrix(NA_real_,nrow(d),2)
for(j in 1:2) residbasis[ii,j]<-resid(lm(bs[ii,j]~d$cmibi_sd[ii]))
d$nl1<-residbasis[,1]; d$nl2<-residbasis[,2]
des<-update(des,nl1=d$nl1[ref],nl2=d$nl2[ref])
sm<-fit('CVD_spline','cvd_death',paste('cmibi_sd + nl1 + nl2 +',m3))
capture.output(regTermTest(sm,~nl1+nl2),file=file.path(out,'spline_nonlinearity.txt'))
capture.output(regTermTest(sm,~cmibi_sd+nl1+nl2),file=file.path(out,'spline_overall.txt'))
for(ev in c('mortstat','cvd_death','cancer_death')) {
  m<-models[[paste0(ev,'_M3')]]
  capture.output(try(cox.zph(m),silent=TRUE),file=file.path(out,paste0(ev,'_PH_diagnostic.txt')))
}
inter<-fit('CVD_sex_interaction','cvd_death',paste('cmibi_sd * factor(riagendr) +',m3))
capture.output(regTermTest(inter,~cmibi_sd:factor(riagendr)),file=file.path(out,'sex_interaction.txt'))
write.csv(do.call(rbind,results),file.path(out,'cox_results.csv'),row.names=FALSE)
counts<-do.call(rbind,lapply(levels(d$cycle),function(s){x<-d[d$cycle==s,]; data.frame(cycle=s,rows=nrow(x),score_complete=sum(is.finite(x$cmibi)),positive_weight=sum(x$weight>0,na.rm=TRUE),eligible=sum(x$eligstat==1,na.rm=TRUE),median_exam_followup=median(x$time,na.rm=TRUE),bp_missing=sum(is.na(x$bp_correct)),dm_missing=sum(is.na(x$dm_correct)),pa_missing=sum(is.na(x$pa_correct)))}))
write.csv(counts,file.path(out,'flow_by_cycle.csv'),row.names=FALSE)
summ<-do.call(rbind,lapply(levels(d$cycle),function(s){x<-d[d$cycle==s&is.finite(d$weight)&d$weight>0,]; data.frame(cycle=s,t(sapply(c('lbxgh','lbxtr','lbdhdd','lbxtc','sbp','dbp','bmxbmi','bmxwaist','crp_mgl','egfr2021','lbxsua','lbxsal','lbxsapsi'),function(v)weighted.mean(x[[v]],x$weight,na.rm=TRUE))))}))
write.csv(summ,file.path(out,'components_by_cycle.csv'),row.names=FALSE)
saveRDS(list(data=d,components=nc,reference=params,models=models),file.path('work','audit_models.rds'))
sink()
