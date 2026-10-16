.libPaths(c('work/Rlib',.libPaths()));Sys.setlocale('LC_CTYPE','Chinese_China.utf8');library(foreign);library(survey);library(survival)
options(survey.lonely.psu='adjust');out<-'outputs/final_analysis';o<-readRDS('work/audit_models.rds');pars<-o$reference;sc<-read.csv('outputs/score_scale.csv')
readx<-function(period,file,cols=NULL){p<-file.path(paste0('E:/NHANES/',period,'年'),paste0(file,'.xpt'));if(!file.exists(p))p<-file.path('work',paste0(file,'.xpt'));z<-read.xport(p);names(z)<-tolower(names(z));if(!is.null(cols))z<-z[,intersect(cols,names(z)),drop=FALSE];stopifnot(!anyDuplicated(z$seqn));z}
pick<-function(x,cols){for(n in cols)if(n%in%names(x))return(x[[n]]);rep(NA_real_,nrow(x))}
yn<-function(x)ifelse(x==1,1,ifelse(x==2,0,NA_real_));frames<-list();flow<-list()
for(i in 1:3){period<-c('1999-2000','2001-2002','2003-2004')[i];s<-c('','_B','_C')[i];file<-function(n)paste0(n,s)
 demo<-readx(period,file('DEMO'),c('seqn','ridageyr','riagendr','ridreth1','dmdeduc2','indfmpir','ridexprg','sdmvstra','sdmvpsu'))
 fs<-list(readx(period,file('BMX'),c('seqn','bmxbmi','bmxwaist')),readx(period,file('BPX'),c('seqn',paste0('bpxsy',1:4),paste0('bpxdi',1:4))),readx(period,if(i==1)'LAB10'else file('L10'),c('seqn','lbxgh')),readx(period,if(i==1)'LAB11'else file('L11'),c('seqn','lbxcrp')),readx(period,if(i==1)'LAB13'else file('L13'),c('seqn','lbxtc','lbdhdl','lbxhdd')),readx(period,if(i==1)'LAB13AM'else file('L13AM'),c('seqn','wtsaf2yr','wtsaf4yr','lbxtr')),readx(period,if(i==1)'LAB18'else file('L40'),c('seqn','lbxsal','lbxsapsi','lbdsapsi','lbxsua','lbxscr','lbdscr')),readx(period,file('SMQ'),c('seqn','smq020','smq040')),readx(period,file('PAQ'),c('seqn','pad200','pad320')),readx(period,file('BPQ'),c('seqn','bpq020','bpq040a','bpq050a','bpq100d')),readx(period,file('DIQ'),c('seqn','diq010','diq050','diq070','did070')))
 d<-demo;for(z in fs)d<-merge(d,z,by='seqn',all.x=TRUE)
 mortdir<-list.dirs('E:/NHANES',recursive=FALSE,full.names=TRUE);mortdir<-mortdir[grepl('1999-2018',mortdir,fixed=TRUE)];lines<-readLines(file.path(mortdir,paste0('NHANES_',gsub('-','_',period),'_MORT_2019_PUBLIC.dat')))
 mort<-data.frame(seqn=as.numeric(substr(lines,1,6)),eligstat=as.numeric(substr(lines,15,15)),mortstat=suppressWarnings(as.numeric(substr(lines,16,16))),cause=suppressWarnings(as.numeric(substr(lines,17,19))),time=suppressWarnings(as.numeric(substr(lines,46,48)))/12)
 d<-merge(d,mort,by='seqn',all.x=TRUE);d$cvd_death<-as.numeric(d$mortstat==1&d$cause%in%c(1,5));d$cancer_death<-as.numeric(d$mortstat==1&d$cause%in%2)
 d<-d[d$ridageyr>=20&!is.na(d$wtsaf2yr)&d$wtsaf2yr>0&(is.na(d$ridexprg)|d$ridexprg!=1),];d$weight<-if(i<3)d$wtsaf4yr*2/3 else d$wtsaf2yr/3
 d$lbdhdd<-pick(d,c('lbdhdl','lbxhdd'));d$lbxsapsi<-pick(d,c('lbxsapsi','lbdsapsi'));d$lbxscr<-pick(d,c('lbxscr','lbdscr'))
 d$sbp<-rowMeans(d[,paste0('bpxsy',1:4)],na.rm=TRUE);d$dbp<-rowMeans(d[,paste0('bpxdi',1:4)],na.rm=TRUE);d$sbp[d$sbp<=50]<-NA;d$dbp[d$dbp<20]<-NA
 scr<-d$lbxscr;if(i==1)scr<-1.013*scr+.147
 k<-ifelse(d$riagendr==1,.9,.7);a<-ifelse(d$riagendr==1,-.302,-.241);d$egfr<-142*pmin(scr/k,1)^a*pmax(scr/k,1)^(-1.2)*.9938^d$ridageyr*ifelse(d$riagendr==2,1.012,1)
 raw<-data.frame(HbA1c=d$lbxgh,TG=log(pmax(d$lbxtr,1)),HDL=-d$lbdhdd,TC=d$lbxtc,SBP=d$sbp,DBP=d$dbp,BMI=d$bmxbmi,waist=d$bmxwaist,CRP=log(pmax(d$lbxcrp*10,.1)),eGFR=-d$egfr,urate=d$lbxsua,albumin=-d$lbxsal,ALP=log(pmax(d$lbxsapsi,1)))
 nc<-Map(function(x,mu,sd)(x-mu)/sd,raw,pars$mean,pars$sd);d$cmibi<-rowSums(as.data.frame(nc));d$cmibi_sd<-(d$cmibi-sc$mean)/sc$sd
 d$age10<-d$ridageyr/10;d$age_square<-d$age10^2;d$race<-factor(d$ridreth1);d$education<-factor(ifelse(d$dmdeduc2%in%1:2,1,ifelse(d$dmdeduc2==3,2,ifelse(d$dmdeduc2%in%4:5,3,NA))))
 d$smoking_cat<-factor(ifelse(d$smq020==2,0,ifelse(d$smq020==1&d$smq040==3,1,ifelse(d$smq020==1&d$smq040%in%1:2,2,NA))),levels=0:2)
 d$pa_correct<-ifelse(d$pad200==1|d$pad320==1,1,ifelse(d$pad200==2&d$pad320==2,0,NA))
 d$bp_correct<-yn(d$bpq050a);skip<-is.na(d$bp_correct)&(d$bpq020==2|d$bpq040a==2);d$bp_correct[which(skip)]<-0
 insulin<-yn(d$diq050);pills<-yn(pick(d,c('did070','diq070')));pills[which(is.na(pills)&d$diq010==2)]<-0;d$dm_correct<-ifelse(insulin==1|pills==1,1,ifelse(insulin==0&pills==0,0,NA));d$lipid_med<-ifelse(d$bpq100d==1,1,0);d$lipid_med[is.na(d$bpq100d)]<-0
 d$period<-period;d$stratum<-factor(d$sdmvstra);d$psu<-interaction(d$sdmvstra,d$sdmvpsu,drop=TRUE)
 flow[[period]]<-data.frame(period=period,adult_positive_fasting_nonpregnant=nrow(d),score_complete=sum(is.finite(d$cmibi_sd)),eligible_score=sum(is.finite(d$cmibi_sd)&d$eligstat==1&d$time>0,na.rm=TRUE))
 frames[[period]]<-d[,c('seqn','period','ridageyr','riagendr','age10','age_square','indfmpir','race','education','smoking_cat','pa_correct','bp_correct','dm_correct','lipid_med','cmibi','cmibi_sd','time','mortstat','cvd_death','cancer_death','eligstat','weight','stratum','psu')]
}
ext<-do.call(rbind,frames);stopifnot(!any(ext$seqn%in%o$data$seqn));ext<-ext[is.finite(ext$cmibi_sd)&ext$eligstat==1&ext$weight>0&ext$time>0,]
write.csv(do.call(rbind,flow),file.path(out,'external_flow.csv'),row.names=FALSE)
des<-svydesign(ids=~psu,strata=~stratum,weights=~weight,data=ext,nest=TRUE);res<-list()
for(ev in c('mortstat','cvd_death','cancer_death'))for(type in c('age_sex','M3')){cov<-if(type=='age_sex')'cmibi_sd+age10+factor(riagendr)'else'cmibi_sd+age10+factor(riagendr)+indfmpir+smoking_cat+education+race+pa_correct+bp_correct+dm_correct+lipid_med';f<-as.formula(paste('Surv(time,',ev,')~',cov));m<-svycoxph(f,design=des,ties='breslow',model=TRUE,x=TRUE,y=TRUE);b<-coef(m)['cmibi_sd'];se<-sqrt(vcov(m)['cmibi_sd','cmibi_sd']);df<-degf(des);res[[paste(ev,type)]]<-data.frame(outcome=ev,model=type,n=m$n,events=m$nevent,HR=exp(b),lower=exp(b-qt(.975,df)*se),upper=exp(b+qt(.975,df)*se),p=2*pt(-abs(b/se),df))}
write.csv(do.call(rbind,res),file.path(out,'external_cox.csv'),row.names=FALSE)
# Locked ten-year prediction models developed only in 2005-2008; all deaths of other causes are observed noncases.
dev<-o$data;dev$age_square<-dev$age10^2;dev<-dev[dev$cycle%in%c('D','E')&dev$eligstat==1&dev$weight>0&dev$time>0&is.finite(dev$cmibi_sd),];dev$y10<-as.numeric(dev$cvd_death==1&dev$time<=10);stopifnot(all(dev$time>=10|dev$mortstat==1));ext$y10<-as.numeric(ext$cvd_death==1&ext$time<=10);stopifnot(all(ext$time>=10|ext$mortstat==1))
f0<-y10~age10+age_square+factor(riagendr);f1<-update(f0,~.+cmibi_sd)
m0<-glm(f0,data=dev,weights=weight/mean(weight),family=quasibinomial());m1<-glm(f1,data=dev,weights=weight/mean(weight),family=quasibinomial());ext$p0<-predict(m0,ext,type='response');ext$p1<-predict(m1,ext,type='response');ext$lp<-qlogis(ext$p1)
auc<-function(y,p,w){ii<-order(p);y<-y[ii];p<-p[ii];w<-w[ii];g<-split(seq_along(p),match(p,unique(p)));cw<-0;num<-0;for(ix in g){w0<-sum(w[ix]*(1-y[ix]));w1<-sum(w[ix]*y[ix]);num<-num+w1*(cw+.5*w0);cw<-cw+w0};num/(sum(w*y)*sum(w*(1-y)))}
metrics<-function(w){ix<-w>0;x<-ext[ix,];x$wb<-w[ix]/mean(w[ix]);c0<-glm(y10~offset(lp),data=x,weights=wb,family=quasibinomial());c1<-glm(y10~lp,data=x,weights=wb,family=quasibinomial());c(AUC_base=auc(x$y10,x$p0,x$wb),AUC_extended=auc(x$y10,x$p1,x$wb),Brier=weighted.mean((x$y10-x$p1)^2,x$wb),calibration_intercept=coef(c0)[1],calibration_slope=coef(c1)[2],observed=weighted.mean(x$y10,x$wb),predicted=weighted.mean(x$p1,x$wb))}
app<-metrics(ext$weight);app<-c(app,AUC_delta=app[2]-app[1]);set.seed(20261015);groups<-split(seq_len(nrow(ext)),ext$stratum,drop=TRUE);boot<-matrix(NA_real_,300,length(app))
for(b in 1:300){mult<-rep(0,nrow(ext));for(ix in groups){ps<-unique(ext$psu[ix]);nh<-length(ps);if(nh<2){mult[ix]<-1;next};picked<-sample(ps,nh-1,replace=TRUE);mult[ix]<-tabulate(match(picked,ps),nbins=nh)[match(ext$psu[ix],ps)]*nh/(nh-1)};mb<-metrics(ext$weight*mult);boot[b,]<-c(mb,AUC_delta=mb[2]-mb[1])}
predres<-data.frame(metric=names(app),estimate=app,lower=apply(boot,2,quantile,.025),upper=apply(boot,2,quantile,.975),development_n=nrow(dev),validation_n=nrow(ext),validation_events=sum(ext$y10));write.csv(predres,file.path(out,'external_locked_prediction.csv'),row.names=FALSE)
write.csv(data.frame(term=names(coef(m1)),coefficient=coef(m1)),file.path(out,'locked_prediction_coefficients.csv'),row.names=FALSE)
ext$decile<-cut(ext$p1,unique(quantile(ext$p1,seq(0,1,.1))),include.lowest=TRUE);cal<-do.call(rbind,lapply(split(ext,ext$decile),function(x)data.frame(n=nrow(x),observed=weighted.mean(x$y10,x$weight),predicted=weighted.mean(x$p1,x$weight))));write.csv(cal,file.path(out,'external_calibration_deciles.csv'),row.names=FALSE)
png(file.path(out,'external_calibration.png'),width=1500,height=1100,res=170);lim<-c(0,max(cal$observed,cal$predicted)*1.1);plot(cal$predicted,cal$observed,xlim=lim,ylim=lim,pch=19,cex=1.2,col='#204963',xlab='Predicted 10-year CVD mortality',ylab='Observed 10-year CVD mortality',main='Locked model in NHANES 1999-2004');abline(0,1,lty=2,col='gray50');lines(cal$predicted,cal$observed,col='#204963');dev.off()
saveRDS(list(external=ext,development=dev,base=m0,extended=m1),'work/external_validation.rds')
