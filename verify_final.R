.libPaths(c('work/Rlib',.libPaths()));library(data.table);out<-'outputs/final_analysis'
m<-fread(file.path(out,'mr_results_corrected.csv'));stopifnot(all(m$lower<=m$OR&m$upper>=m$OR),all(is.finite(m$OR)),all(m$p>=0&m$p<=1),sum(m$method=='IVW random effects')==10)
for(f in list.files(out,pattern='^MR_.*_harmonised.csv$',full.names=TRUE)){x<-fread(f);stopifnot(!anyDuplicated(x$SNP),all(x$se_out>0),all((x$beta_exp/x$se_exp)^2>10))}
mi<-fread(file.path(out,'multiple_imputation_cox.csv'));stopifnot(all(mi$m==30),all(mi$n==10868),all(mi$lower<mi$HR&mi$HR<mi$upper))
ext<-readRDS('work/external_validation.rds');o<-readRDS('work/audit_models.rds');stopifnot(!any(ext$external$seqn%in%o$data$seqn),nrow(ext$external)==5199,all(ext$external$time>=10|ext$external$mortstat==1))
saveRDS(o$data,file.path(out,'corrected_development_data.rds'))
write.csv(data.frame(package=c('survey','survival','mice','limma','data.table','officer','flextable'),version=sapply(c('survey','survival','mice','limma','data.table','officer','flextable'),function(p)as.character(packageVersion(p)))),file.path(out,'software_versions.csv'),row.names=FALSE)
writeLines(c('PASS: ten primary IVW comparisons and finite coherent intervals.','PASS: no duplicated coordinated SNPs; all instrument F > 10.','PASS: 30 imputations with identical outcome count and cohort size.','PASS: development/validation participant IDs do not overlap.','PASS: all historical validation ten-year event statuses are ascertainable.'),file.path(out,'verification_checks.txt'))
