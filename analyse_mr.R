.libPaths(c('work/Rlib',.libPaths()));library(data.table);source('work/scripts/mr_methods.R');source('work/add_rsq.R');source('work/mr_presso.R')
out<-'outputs/final_analysis';set.seed(20261015)
# Validate the optimized algebra against original MR-PRESSO on identical small simulated inputs.
v<-data.frame(bx=seq(.05,.15,length.out=8),by=seq(.05,.15,length.out=8)*.3+rnorm(8,0,.015),sx=rep(.01,8),sy=rep(.02,8))
orig<-suppressWarnings(mr_presso('by','bx','sy','sx',data=v,OUTLIERtest=TRUE,DISTORTIONtest=FALSE,NbDistribution=200,seed=765))
fast<-presso_fast(v$bx,v$by,v$sx,v$sy,B=200,seed=765)
stopifnot(abs(orig[['MR-PRESSO results']][['Global Test']]$RSSobs-fast$RSS)<1e-8,abs(orig[['MR-PRESSO results']][['Global Test']]$Pvalue-fast$global_p)<1e-10)
writeLines('Optimized univariable residual algorithm matched original MR-PRESSO observed RSS and Monte Carlo global P on a fixed-seed 200-replicate validation example.',file.path(out,'mr_presso_validation.txt'))
harm<-function(x,y){m<-merge(x,y,by='SNP',suffixes=c('_exp','_out'));m<-m[!duplicated(SNP)];same<-m$EA_exp==m$EA_out&m$OA_exp==m$OA_out;flip<-m$EA_exp==m$OA_out&m$OA_exp==m$EA_out
 complement<-function(s)chartr('ACGT','TGCA',s);cs<-m$EA_exp==complement(m$EA_out)&m$OA_exp==complement(m$OA_out);cf<-m$EA_exp==complement(m$OA_out)&m$OA_exp==complement(m$EA_out)
 keep<-same|flip|cs|cf;m<-m[keep];sg<-ifelse((flip|cf)[keep],-1,1);m$beta_out<-m$beta_out*sg;m$EAF_out<-ifelse(sg==-1,1-m$EAF_out,m$EAF_out);m<-m[is.finite(beta_exp)&is.finite(beta_out)&is.finite(se_exp)&is.finite(se_out)&se_exp>0&se_out>0];m}
pairs<-rbind(data.frame(exposure=rep(c('TG','HDL'),each=2),outcome=rep(c('CHD','HF'),2)),data.frame(exposure=c('BMI','CRP','SBP','DBP','eGFR','urate'),outcome='stroke'))
res<-list();diag<-list();steg<-list();looall<-list()
for(j in 1:nrow(pairs)){ex<-pairs$exposure[j];oc<-pairs$outcome[j];x<-fread(file.path('work/mr',paste0(ex,'_instruments.csv')));y<-fread(file.path('work/mr',paste0(oc,'_lookup.csv')));m<-harm(x,y);cat(ex,oc,nrow(m),'\n');flush.console();bx<-m$beta_exp;by<-m$beta_out;sx<-m$se_exp;sy<-m$se_out;n<-nrow(m)
 if(n<4)next
 est<-ivw(bx,by,sy);med<-wmedian(by/bx,bx^2/sy^2);set.seed(20261015);mb<-replicate(2000,{xx<-rnorm(n,bx,sx);yy<-rnorm(n,by,sy);wmedian(yy/xx,xx^2/sy^2)});ms<-sd(mb)
 sg<-sign(bx);eg<-lm(by*sg~I(abs(bx)),weights=1/sy^2);es<-summary(eg)$coef;eb<-es[2,1];ese<-es[2,2];ep<-es[2,4]
 make<-function(method,b,se,p)data.frame(exposure=ex,outcome=oc,method=method,n_snp=n,beta=b,se=se,OR=exp(b),lower=exp(b-1.96*se),upper=exp(b+1.96*se),p=p)
 res[[length(res)+1]]<-make('IVW random effects',est['beta'],est['se'],est['p']);res[[length(res)+1]]<-make('Weighted median bootstrap',med,ms,2*pnorm(-abs(med/ms)));res[[length(res)+1]]<-make('MR Egger',eb,ese,ep)
 # Steiger comparison uses sample-size-derived quantitative correlations and log-odds-derived outcome correlations.
 nn<-m$N_exp;if(ex=='CRP')nn[!is.finite(nn)]<-575531
 rx2<-(bx/sx)^2/((bx/sx)^2+nn-2)
 nc<-if(oc=='CHD')56650 else if(oc=='HF')37653 else 73652
 nt<-if(oc%in%c('CHD','HF'))500348 else 1308460
 prev<-if(oc=='CHD').05 else if(oc=='HF').02 else .03
 ry<-get_r_from_lor(by,m$EAF_out,ncase=nc,ncontrol=nt-nc,prevalence=prev)
 keep<-rx2>ry^2&is.finite(rx2)&is.finite(ry);m$steiger_keep<-keep;m$R2_exposure<-rx2;m$R2_outcome<-ry^2
 # Population prevalence sensitivity: 1%, 5%, 10%; direction counts, not causal proof.
 for(pr in c(.01,.05,.1)){r2<-get_r_from_lor(by,m$EAF_out,nc,nt-nc,pr)^2;steg[[length(steg)+1]]<-data.frame(exposure=ex,outcome=oc,prevalence=pr,n=n,kept=sum(rx2>r2,na.rm=TRUE))}
 if(sum(keep)>=4){z<-ivw(bx[keep],by[keep],sy[keep]);r<-make('IVW after Steiger',z['beta'],z['se'],z['p']);r$n_snp<-sum(keep);res[[length(res)+1]]<-r}
 pz<-presso_fast(bx,by,sx,sy,B=20000,seed=20261015);m$presso_outlier_p<-pz$outlier_p
 if(all(is.finite(pz$corrected[1:3]))){r<-make('MR PRESSO outlier corrected',pz$corrected[1],pz$corrected[2],pz$corrected[3]);r$n_snp<-n-length(pz$outliers);res[[length(res)+1]]<-r}
 diag[[length(diag)+1]]<-data.frame(exposure=ex,outcome=oc,n_snp=n,F_min=min((bx/sx)^2),F_mean=mean((bx/sx)^2),Q=est['Q'],Q_p=est['Q_p'],egger_intercept=es[1,1],egger_intercept_p=es[1,4],steiger_kept=sum(keep),PRESSO_global_p=pz$global_p,PRESSO_outliers=length(pz$outliers),PRESSO_distortion_pct=100*pz$distortion,PRESSO_distortion_p=pz$distortion_p,simulations=20000)
 fwrite(m,file.path(out,paste0('MR_',ex,'_',oc,'_harmonised.csv')))
 lo<-do.call(rbind,lapply(seq_len(n),function(k){z<-ivw(bx[-k],by[-k],sy[-k]);data.frame(exposure=ex,outcome=oc,SNP=m$SNP[k],OR=exp(z['beta']),lower=exp(z['beta']-1.96*z['se']),upper=exp(z['beta']+1.96*z['se']))}));looall[[length(looall)+1]]<-lo
 png(file.path(out,paste0('MR_',ex,'_',oc,'_scatter.png')),width=1400,height=950,res=170);plot(bx,by,pch=19,col='#4d7187',xlab=paste(ex,'GWAS effect'),ylab=paste(oc,'log odds effect'),main=paste(ex,'to',oc));abline(0,est['beta'],col='#a44648',lwd=2);abline(es[1,1],eb,col='#9a873f',lty=2);legend('topleft',c('IVW','MR Egger'),col=c('#a44648','#9a873f'),lty=c(1,2),bty='n');dev.off()
}
r<-rbindlist(res);r$FDR_IVW<-NA_real_;ii<-r$method=='IVW random effects';r$FDR_IVW[ii]<-p.adjust(r$p[ii],'BH');fwrite(r,file.path(out,'mr_results_corrected.csv'));fwrite(rbindlist(diag),file.path(out,'mr_diagnostics.csv'));fwrite(rbindlist(steg),file.path(out,'steiger_prevalence_sensitivity.csv'));fwrite(rbindlist(looall),file.path(out,'mr_leave_one_out.csv'))
z<-r[method=='IVW random effects'];png(file.path(out,'mr_forest.png'),width=1800,height=1300,res=180);par(mar=c(5,12,3,2));yy<-rev(seq_len(nrow(z)));plot(z$OR,yy,log='x',xlim=range(c(z$lower,z$upper)),yaxt='n',pch=19,xlab='Outcome odds ratio per published GWAS exposure unit',ylab='',main='LD-clumped two-sample MR');segments(z$lower,yy,z$upper,yy,lwd=2);axis(2,at=yy,labels=paste(z$exposure,'to',z$outcome),las=1);abline(v=1,lty=2,col='gray50');dev.off()
saveRDS(list(results=r,diagnostics=rbindlist(diag)),'work/mr_results.rds')
