

aa <- read.table("Hutterite_taxa_abundance_table_genus_level2.txt")
aa[37,1] <- "Lachnospiraceae"
aa1 <- log(aa[,-c(1:5)]+1)

aa1.n <- colnames(aa1)
s.i <- c()
for(i in 1:length(aa1.n)){
  
  tmp1 <- strsplit(aa1.n[i],"_")[[1]][3]
  tmp2 <- strsplit(tmp1,"")[[1]][1]
  s.i <- c(s.i,tmp2)
}

sum.i <- which(s.i=="S")
win.i <- which(s.i=="W")

gna <- rownames(aa)

gna.1 <- c()
for(i in 1:length(gna)){
  
  gna.1 <- c(gna.1,strsplit(gna[i],"_")[[1]][2])
}

gs <- aa1[,sum.i]
gw <- aa1[,win.i]

gs.sum <- rowSums(gs)
gw.sum <- rowSums(gw)


res <- list(gall=aa1,gs=gs,gw=gw,gs.sum=gs.sum,gw.sum=gw.sum,gname=gna)

save(res,file="dat.RData")
