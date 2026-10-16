.libPaths(c('work/Rlib',.libPaths()));library(data.table)
for(f in c('format_function.R','strength_mvmr.R','phenocov_mvmr.R'))source(file.path('work/MVMR_source/MVMR-master/R',f))
out<-'outputs/final_analysis';ids<-readLines('work/mr/lipid_union_snps.txt');ref<-fread('work/mr/HDL_lookup.csv');ref<-ref[SNP%in%ids][!duplicated(SNP)];ref$EA0<-ref$EA;ref$OA0<-ref$OA
nm<-c('HDL','TG','LDL');B<-S<-matrix(NA_real_,nrow(ref),3,dimnames=list(NULL,nm))
align<-function(r,x){i<-match(r$SNP,x$SNP);x<-x[i];sg<-ifelse(r$EA0==x$EA&r$OA0==x$OA,1,ifelse(r$EA0==x$OA&r$OA0==x$EA,-1,NA_real_));list(beta=x$beta*sg,se=x$se)}
for(k in 1:3){x<-fread(file.path('work/mr',paste0(nm[k],'_lookup.csv')));a<-align(ref,x);B[,k]<-a$beta;S[,k]<-a$se}
res<-list();strength<-list()
for(oc in c('CHD','HF')){y<-align(ref,fread(file.path('work/mr',paste0(oc,'_lookup.csv'))));ok<-complete.cases(B,S,y$beta,y$se);X<-B[ok,,drop=FALSE];seX<-S[ok,,drop=FALSE];by<-y$beta[ok];sy<-y$se[ok];w<-1/sy^2
 bh<-solve(crossprod(X,X*w),crossprod(X,by*w));Q<-sum(w*(by-X%*%bh)^2);V<-solve(crossprod(X,X*w))*max(1,Q/(nrow(X)-3));se<-sqrt(diag(V));res[[oc]]<-data.frame(outcome=oc,exposure=nm,n_snp=nrow(X),beta=as.numeric(bh),se=se,OR=exp(bh),lower=exp(bh-1.96*se),upper=exp(bh+1.96*se),p=2*pnorm(-abs(bh/se)),Q=Q,Q_p=pchisq(Q,nrow(X)-3,lower.tail=FALSE))
 inp<-format_mvmr(X,by,seX,sy,RSID=ref$SNP[ok]);fs<-suppressWarnings(strength_mvmr(inp,gencov=0));strength[[paste0(oc,'_zero')]]<-data.frame(outcome=oc,covariance='zero',exposure=nm,F_conditional=as.numeric(fs[1,]))
 # External phenotypic correlations estimated from NHANES non-Hispanic White fasting adults; covariance sensitivity only.
 o<-readRDS('work/audit_models.rds');d<-o$data;z<-data.frame(HDL=d$lbdhdd,TG=log(d$lbxtr),LDL=d$lbxtc-d$lbdhdd-d$lbxtr/5);cc<-d$ridreth1==3&d$lbxtr<400&d$weight>0&complete.cases(z);rho<-cor(z[cc,]);gc<-phenocov_mvmr(rho,seX);fp<-strength_mvmr(inp,gencov=gc);strength[[paste0(oc,'_rho')]]<-data.frame(outcome=oc,covariance='NHANES phenotypic correlation sensitivity',exposure=nm,F_conditional=as.numeric(fp[1,]));write.csv(rho,file.path(out,'lipid_correlation_sensitivity.csv'))
 fwrite(data.table(SNP=ref$SNP[ok],X,seX,beta_out=by,se_out=sy),file.path(out,paste0('MVMR_',oc,'_harmonised.csv')))
}
fwrite(rbindlist(res),file.path(out,'multivariable_mr.csv'));fwrite(rbindlist(strength),file.path(out,'multivariable_conditional_F.csv'))
