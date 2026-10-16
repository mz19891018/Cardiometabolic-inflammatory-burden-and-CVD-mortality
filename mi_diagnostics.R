.libPaths(c('work/Rlib',.libPaths()));library(mice)
m<-readRDS('work/mice_object.rds')
# Conventional between-chain/within-chain R-hat on iteration 11-20 chain means.
cm<-m$chainMean;res<-list()
for(v in dimnames(cm)[[1]]){a<-cm[v,11:20,,drop=FALSE];a<-matrix(a,nrow=10,ncol=m$m);if(!all(is.finite(a)))next;W<-mean(apply(a,2,var));B<-10*var(colMeans(a));rh<-sqrt(((9/10)*W+B/10)/W);res[[v]]<-data.frame(variable=v,Rhat_chain_means=rh,iterations=10,chains=m$m)}
write.csv(do.call(rbind,res),'outputs/final_analysis/mi_convergence_statistics.csv',row.names=FALSE)
capture.output(print(do.call(rbind,res)),file='outputs/final_analysis/mi_final_convergence.txt');ev<-m$loggedEvents;write.csv(as.data.frame(table(ev$out)),'outputs/final_analysis/mi_logged_event_summary.csv',row.names=FALSE)
png('outputs/final_analysis/mi_convergence.png',width=1800,height=1600,res=170);par(mfrow=c(4,2),mar=c(3,4,2,1))
for(v in c('indfmpir','smoking_cat','education','pa_correct')){a<-matrix(m$chainMean[v,,],nrow=20);b<-sqrt(matrix(m$chainVar[v,,],nrow=20));matplot(1:20,a,type='l',lty=1,col=rainbow(m$m),xlab='Iteration',ylab='Imputed mean',main=v);matplot(1:20,b,type='l',lty=1,col=rainbow(m$m),xlab='Iteration',ylab='Imputed SD',main=v)};dev.off()
