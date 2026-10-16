.libPaths(c('work/Rlib',.libPaths())); library(mice); library(survey); library(survival)
options(survey.lonely.psu='adjust'); set.seed(20261015)
out<-'outputs/final_analysis'; dir.create(out,recursive=TRUE,showWarnings=FALSE)
o<-readRDS('work/audit_models.rds'); d<-o$data
d<-d[d$eligstat==1&d$weight>0&is.finite(d$cmibi_sd)&d$time>0,]
haz<-function(ev){fit<-survfit(Surv(d$time,d[[ev]])~1); h<-cumsum(fit$n.event/fit$n.risk); c(0,h)[findInterval(d$time,fit$time)+1]}
d$H_all<-haz('mortstat');d$H_cvd<-haz('cvd_death');d$H_cancer<-haz('cancer_death')
v<-c('age10','riagendr','indfmpir','smoking_cat','education','race','pa_correct','bp_correct','dm_correct','lipid_med','cmibi_sd','time','mortstat','cvd_death','cancer_death','H_all','H_cvd','H_cancer','cycle')
x<-d[,v]; x$riagendr<-factor(x$riagendr); for(n in c('pa_correct','bp_correct','dm_correct'))x[[n]]<-factor(x[[n]],levels=0:1)
x$logweight<-log(d$weight); x$age_square<-x$age10^2
ini<-mice(x,maxit=0,printFlag=FALSE); method<-ini$method; method['indfmpir']<-'pmm'; method[c('smoking_cat','education')]<-'polyreg'; method[c('pa_correct','bp_correct','dm_correct')]<-'logreg'
pred<-ini$predictorMatrix; pred[,c('time')]<-0; pred[c('H_all','H_cvd','H_cancer','age_square','logweight'),]<-0
imp<-mice(x,m=30,maxit=20,method=method,predictorMatrix=pred,seed=20261015,printFlag=FALSE)
saveRDS(imp,'work/mice_object.rds'); if(!is.null(imp$loggedEvents))write.csv(imp$loggedEvents,file.path(out,'mi_logged_events.csv'),row.names=FALSE)
write.csv(data.frame(variable=names(x),missing=colSums(is.na(x)),pct=colMeans(is.na(x))*100),file.path(out,'missingness.csv'),row.names=FALSE)
res<-list(); each<-list()
for(ev in c('mortstat','cvd_death','cancer_death')) {
 b<-u<-numeric(30)
 for(i in 1:30){z<-complete(imp,i);for(n in names(z))d[[n]]<-z[[n]];des<-svydesign(ids=~psu,strata=~stratum,weights=~weight,data=d,nest=TRUE)
 f<-as.formula(paste('Surv(time,',ev,')~cmibi_sd+age10+factor(riagendr)+indfmpir+smoking_cat+education+race+pa_correct+bp_correct+dm_correct+lipid_med'))
 fit<-svycoxph(f,design=des,ties='breslow',model=TRUE,x=TRUE,y=TRUE);b[i]<-coef(fit)['cmibi_sd'];u[i]<-vcov(fit)['cmibi_sd','cmibi_sd'] }
 pooled<-pool.scalar(b,u,n=78,k=1);se<-sqrt(pooled$t);df<-pooled$df;bt<-pooled$qbar
 res[[ev]]<-data.frame(outcome=ev,n=nrow(d),events=sum(d[[ev]]),m=30,HR=exp(bt),lower=exp(bt-qt(.975,df)*se),upper=exp(bt+qt(.975,df)*se),p=2*pt(-abs(bt/se),df),df=df,FMI=pooled$fmi,between=pooled$b,within=pooled$ubar)
 each[[ev]]<-data.frame(outcome=ev,imputation=1:30,beta=b,variance=u)
}
write.csv(do.call(rbind,res),file.path(out,'multiple_imputation_cox.csv'),row.names=FALSE);write.csv(do.call(rbind,each),file.path(out,'multiple_imputation_estimates.csv'),row.names=FALSE)
png(file.path(out,'mi_convergence.png'),width=1700,height=1100,res=170);print(plot(imp,c('indfmpir','smoking_cat','education','pa_correct')));dev.off()
capture.output(imp,file=file.path(out,'mi_diagnostics.txt'))
