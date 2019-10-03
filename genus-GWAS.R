

load("dat.RData")
load("xuhao.RData")
load("snp_3_1.RData")
source("genus-sup.R")
load("Pheno.RData")


require(foreign);
require(nnet);
library(MASS);
library(lmtest)
library(parallel)

pheno.mut <- all.net(res,method="mut",xuhao)
pheno.ant <- all.net(res,method="ant",xuhao)
pheno.agg <- all.net(res,method="agg",xuhao)
pheno.alt <- all.net(res,method="alt",xuhao)

pheno <- rbind(pheno.mut,pheno.ant,pheno.agg,pheno.alt)

save(pheno,file="Pheno.RData")

pheno.name <- c(paste("mut_",rownames(pheno)[1:12],sep=""),paste("ant_",rownames(pheno)[1:12],sep=""),
                paste("agg_",rownames(pheno)[1:12],sep=""),paste("alt_",rownames(pheno)[1:12],sep=""))

for(i in 1:dim(pheno)[1]){
  ret <-(gwas.scan(snp=nsnp[,1:88],trait=pheno[i,],n.cores=80,proc=pv_continuous))
  file.name <- paste(pheno.name[i],".RData",sep="")
  save(ret,file=file.name)
  cat("i=",i,pheno.name[i],"\n")
}


pv.t <- c()
for(i in 1:100){
  np <- sample(1:88,88)
  ret1 <- as.numeric(gwas.scan(snp=nsnp[np,1:88],trait=pheno.mut$Con.s,n.cores=24,proc=pv_continuous))
  pv.t <- c(pv.t,max(-log(ret1),na.rm=T))
  cat("i=",i,"\n")
}
