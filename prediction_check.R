.libPaths(c('work/Rlib',.libPaths()))
library(survey); library(survival)
obj<-readRDS('work/audit_models.rds'); d<-obj$data
set.seed(20261015)
# Ten-year absolute CVD mortality: only 2005-2008 recruitment, with full follow-up support.
# Other-cause death before year 10 is an observed non-case, never a censoring exclusion.
d$y10<-as.numeric(d$cvd_death==1 & d$time<=10)
d$known10<-d$time>=10|d$mortstat==1
d$mets_correct<-as.numeric((d$bmxwaist>=ifelse(d$riagendr==1,102,88))) + as.numeric(d$lbxtr>=150) + as.numeric(d$lbdhdd<ifelse(d$riagendr==1,40,50)) + as.numeric(d$sbp>=130|d$dbp>=85|d$bp_correct==1) + as.numeric(d$lbxglu>=100|d$dm_correct==1)
# AIP uses mmol/L; a constant shift does not change a fitted continuous regression effect.
d$aip_molar<-log10((d$lbxtr/88.57)/(d$lbdhdd/38.67))
vars<-c('y10','age10','riagendr','indfmpir','smoking_cat','education','race','pa_correct','bp_correct','dm_correct','lipid_med','cmibi_sd','mets_correct','aip_molar','tyg','weight')
dd<-d[d$cycle%in%c('D','E')&d$eligstat==1&d$time>0&d$weight>0&d$known10&complete.cases(d[,vars]),]
stopifnot(all(dd$known10),all(is.finite(dd$weight)))
base<-~age10+factor(riagendr)+indfmpir+smoking_cat+education+race+pa_correct+bp_correct+dm_correct+lipid_med
X<-model.matrix(base,dd); y<-dd$y10; w<-dd$weight/mean(dd$weight)
markers<-c(CMIBI='cmibi_sd',MetS='mets_correct',AIP='aip_molar',TyG='tyg')
metric<-function(a,b,y,w) {
  em<-function(x,g)sum(x[g]*w[g])/sum(w[g])
  sig<-sign(b-a); ev<-y==1; ne<-y==0
  c(NRI=em(sig,ev)-em(sig,ne),IDI=(em(b,ev)-em(b,ne))-(em(a,ev)-em(a,ne)))
}
pred<-function(x,w){ f<-suppressWarnings(glm.fit(x,y,weights=w,family=quasibinomial())); if(!f$converged)stop('GLM convergence failure'); as.numeric(plogis(x%*%f$coefficients)) }
estimate<-function(wb) {a<-pred(X,wb); sapply(markers,function(v)metric(a,pred(cbind(X,dd[[v]]),wb),y,wb))}
app<-estimate(w)
boots<-array(NA_real_,c(2,4,300)); groups<-split(seq_len(nrow(dd)),dd$stratum,drop=TRUE)
for(b in 1:300) {
  mult<-rep(0,nrow(dd))
  for(ix in groups) { ps<-unique(dd$psu[ix]); nh<-length(ps); if(nh<2)stop('Single PSU stratum'); picked<-sample(ps,nh-1,replace=TRUE); mult[ix]<-tabulate(match(picked,ps),nbins=nh)[match(dd$psu[ix],ps)]*nh/(nh-1) }
  boots[,,b]<-tryCatch(estimate(w*mult),error=function(e)matrix(NA_real_,2,4))
}
res<-do.call(rbind,lapply(seq_along(markers),function(j)data.frame(marker=names(markers)[j],n=nrow(dd),events=sum(y),NRI=app[1,j],NRI_low=quantile(boots[1,j,],.025,na.rm=TRUE),NRI_high=quantile(boots[1,j,],.975,na.rm=TRUE),IDI=app[2,j],IDI_low=quantile(boots[2,j,],.025,na.rm=TRUE),IDI_high=quantile(boots[2,j,],.975,na.rm=TRUE),valid_boot=sum(is.finite(boots[1,j,])))))
write.csv(res,'outputs/nri_early_cycles_apparent.csv',row.names=FALSE)
# This is apparent, internally fitted discrimination. It is not external validation or optimism correction.
capture.output(print(res),file='outputs/prediction_notes.txt')
