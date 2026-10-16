.libPaths(c('work/Rlib',.libPaths())); library(survival)
o<-readRDS('work/audit_models.rds'); d<-o$data
cv<-as.matrix(d[,paste0('mcq160',letters[2:6])]); no_cvd<-rowSums(cv==2,na.rm=TRUE)==5
vars<-c('time','cvd_death','age10','riagendr','sbp','lbxtc','lbdhdd','diabetes_correct','smoking','egfr2021','bp_correct','lipid_med','cmibi_sd','weight')
dd<-droplevels(d[no_cvd&d$ridageyr>=30&d$ridageyr<=79&d$eligstat==1&d$time>0&d$weight>0&complete.cases(d[,vars]),])
dd$smoker<-as.numeric(dd$smoking==2)
f0<-Surv(time,cvd_death)~age10+factor(riagendr)+sbp+lbxtc+lbdhdd+diabetes_correct+smoker+egfr2021+bp_correct+lipid_med
f1<-update(f0,~.+cmibi_sd)
calc<-function(w) {
 x<-dd[w>0,]; x$wb<-w[w>0]/mean(w[w>0]); m0<-coxph(f0,data=x,weights=wb,ties='breslow'); m1<-coxph(f1,data=x,weights=wb,ties='breslow')
 x$r0<-predict(m0,type='lp'); x$r1<-predict(m1,type='lp')
 c0<-concordance(Surv(time,cvd_death)~r0,data=x,weights=wb,reverse=TRUE,timewt='n')$concordance
 c1<-concordance(Surv(time,cvd_death)~r1,data=x,weights=wb,reverse=TRUE,timewt='n')$concordance
 c(base=c0,extended=c1,delta=c1-c0)
}
app<-calc(dd$weight); set.seed(20261015); groups<-split(seq_len(nrow(dd)),dd$stratum,drop=TRUE); boot<-matrix(NA_real_,300,3)
for(b in 1:300){mult<-rep(0,nrow(dd));for(ix in groups){ps<-unique(dd$psu[ix]);nh<-length(ps);if(nh<2)stop('single PSU');picked<-sample(ps,nh-1,replace=TRUE);mult[ix]<-tabulate(match(picked,ps),nbins=nh)[match(dd$psu[ix],ps)]*nh/(nh-1)};boot[b,]<-tryCatch(calc(dd$weight*mult),error=function(e)rep(NA_real_,3))}
write.csv(data.frame(n=nrow(dd),events=sum(dd$cvd_death),C_base=app[1],C_extended=app[2],delta=app[3],delta_lower=quantile(boot[,3],.025,na.rm=TRUE),delta_upper=quantile(boot[,3],.975,na.rm=TRUE),valid_boot=sum(is.finite(boot[,3]))),'outputs/discrimination_apparent.csv',row.names=FALSE)
