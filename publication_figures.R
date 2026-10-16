.libPaths(c('work/Rlib',.libPaths()));library(ggplot2);library(ggrepel)
out<-'outputs/final_analysis'
for(nm in c('GSE98895_MetS','GSE9006_T1D','GSE9006_T2D')){d<-read.csv(file.path(out,paste0(nm,'_limma.csv')));d$mlog<- -log10(pmax(d$adj.P.Val,1e-300));d$sig<-d$adj.P.Val<.05&abs(d$logFC)>=.5;lab<-head(d[d$sig,],8)
 p<-ggplot(d,aes(logFC,mlog))+geom_point(aes(color=sig),size=.65,alpha=.8)+scale_color_manual(values=c('#8c9caa','#b13f4c'),guide='none')+geom_vline(xintercept=c(-.5,.5),linetype=2,color='gray65',linewidth=.3)+geom_hline(yintercept=-log10(.05),linetype=2,color='gray65',linewidth=.3)+geom_text_repel(data=lab,aes(label=gene),seed=20261015,size=3,box.padding=.5,point.padding=.3,max.overlaps=Inf,min.segment.length=0,color='black')+labs(title=paste(nm,'limma'),x='Adjusted log2 fold change',y='-log10 FDR')+scale_y_continuous(expand=expansion(mult=c(.02,.15)))+theme_classic(base_size=11)+theme(plot.title=element_text(face='bold',hjust=.5),plot.margin=margin(12,18,12,12))
 ggsave(file.path(out,paste0(nm,'_volcano.png')),p,width=9,height=6.3,dpi=180,bg='white')
}
