ivw<-function(bx,by,sy){w<-1/sy^2;b<-sum(w*bx*by)/sum(w*bx^2);Q<-sum(w*(by-b*bx)^2);se<-sqrt(max(1,Q/(length(bx)-1))/sum(w*bx^2));c(beta=b,se=se,p=2*pnorm(-abs(b/se)),Q=Q,Q_p=pchisq(Q,length(bx)-1,lower.tail=FALSE))}
wmedian<-function(x,w){j<-order(x);x<-x[j];w<-w[j]/sum(w);p<-cumsum(w)-w/2;approx(p,x,xout=.5,rule=2)$y}
# Algebraic single-exposure implementation of the public MR-PRESSO residual algorithm.
# Uses the same leave-one-out slopes, exposure/outcome simulations, RSS, and Bonferroni test.
presso_fast<-function(bx,by,sx,sy,B=20000,seed=20261015){set.seed(seed);sg<-sign(bx);bx<-bx*sg;by<-by*sg;n<-length(bx);w<-1/sy^2
 loo<-function(x,y){(sum(w*x*y)-w*x*y)/(sum(w*x*x)-w*x*x)}
 b0<-loo(bx,by);dif<-by-bx*b0;rss<-sum(w*dif^2);rssp<-numeric(B);outcnt<-numeric(n)
 for(i in 1:B){x<-rnorm(n,bx,sx);y<-rnorm(n,bx*b0,sy);br<-loo(x,y);rssp[i]<-sum(w*(y-x*br)^2);outcnt<-outcnt+as.numeric((y-x*b0)^2>dif^2)}
 global<-mean(rssp>rss);op<-pmin(1,(outcnt/B)*n);outs<-if(global<.05)which(op<=.05)else integer();raw<-ivw(bx,by,sy);cor<-if(length(outs)>0&&length(outs)<n-3)ivw(bx[-outs],by[-outs],sy[-outs])else rep(NA_real_,5)
 distortion<-NA_real_;dp<-NA_real_
 if(length(outs)>0&&length(outs)<n-3){distortion<-(raw['beta']-cor['beta'])/abs(cor['beta']);bias<-numeric(B);for(j in 1:B){ii<-c(outs,sample(setdiff(seq_len(n),outs),n-length(outs),replace=TRUE));ii<-ii[seq_len(n-length(outs))];b<-sum(w[ii]*bx[ii]*by[ii])/sum(w[ii]*bx[ii]^2);bias[j]<-(raw['beta']-b)/abs(b)};dp<-mean(abs(bias)>abs(distortion))}
 list(global_p=global,outlier_p=op,outliers=outs,raw=raw,corrected=cor,distortion=distortion,distortion_p=dp,RSS=rss,B=B)
}
