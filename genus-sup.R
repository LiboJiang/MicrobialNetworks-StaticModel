
all.net.geno <- function(res,xuhao,geno){
  
  ds <- res$gs
  dw <- res$gw
  
  sname <- colnames(ds)
  sname1 <- c()
  for(i in 1:length(sname)){
    sname1 <- c(sname1,strsplit(sname[i],"_S")[[1]][1])
  }
  sname2 <- unique(sname1)
  
  
  wname <- colnames(dw)
  wname1 <- c()
  for(i in 1:length(wname)){
    wname1 <- c(wname1,strsplit(wname[i],"_W")[[1]][1])
  }
  wname2 <- unique(wname1)
  
  sw <- table(c(sname2,wname2))
  swname <- names(sw)
  swname1 <- swname[which(sw==2)]
  
  ds1 <- c()
  for(i in 1:length(swname1)){
    
    index <- which(sname1==swname1[i])
    ds1 <- cbind(ds1,as.numeric(rowSums(ds[,index])))
  }
  colnames(ds1) <- swname1
  rownames(ds1) <- res$gname
  
  dw1 <- c()
  for(i in 1:length(swname1)){
    
    index <- which(wname1==swname1[i])
    dw1 <- cbind(dw1,as.numeric(rowSums(dw[,index])))
  }
  colnames(dw1) <- swname1
  rownames(dw1) <- res$gname
  
  ds1 <- ds1[,xuhao]
  dw1 <- dw1[,xuhao]
  nnn <- length(swname1)
  AA.i <- which(geno==2)
  Aa.i <- which(geno==1)
  aa.i <- which(geno==0)
  
  s.AAn <- thr.compute(pheno=rowSums(ds1[,AA.i]))
  s.Aan <- thr.compute(pheno=rowSums(ds1[,Aa.i]))
  s.aan <- thr.compute(pheno=rowSums(ds1[,aa.i]))
  w.AAn <- thr.compute(pheno=rowSums(dw1[,AA.i]))
  w.Aan <- thr.compute(pheno=rowSums(dw1[,Aa.i]))
  w.aan <- thr.compute(pheno=rowSums(dw1[,aa.i]))
  
  astmp_AA <- network.sparse(gg=s.AAn,en=200)
  astmp_Aa <- network.sparse(gg=s.Aan,en=200)
  astmp_aa <- network.sparse(gg=s.aan,en=200)
  
  awtmp_AA <- network.sparse(gg=w.AAn,en=200)
  awtmp_Aa <- network.sparse(gg=w.Aan,en=200)
  awtmp_aa <- network.sparse(gg=w.aan,en=200)
  
  alist <- list(astmp_AA=astmp_AA,astmp_Aa=astmp_Aa,astmp_aa=astmp_aa,
                awtmp_AA=awtmp_AA,awtmp_Aa=awtmp_Aa,awtmp_aa=awtmp_aa)
  return (alist)
}


gwas.scan <- function(snp,trait,n.cores,proc){
  
  nsnp <- dim(snp)[1]
  
  grp <- floor(nsnp/n.cores)
  grp.i <- c()
  if(n.cores==1){
    grp.i <- c(grp.i,rep(1,nsnp))
  }else{
    for(ii in 1:n.cores){
      if(ii==n.cores){
        grp.i <- c(grp.i,rep(ii,nsnp-grp*(ii-1)))
      }else{
        grp.i <- c(grp.i,rep(ii,grp))
      }
    }
  }
  
  grp.ii <- unique(grp.i)
  
  res.list <- mclapply(grp.ii, function(i)
  {
    snp.c <- 	which(grp.i==i)
    A <- sapply(snp.c, proc, snp=snp, trait=trait);
    return (A);
  }, mc.cores=n.cores )
  
  res1 <- do.call("c", res.list)
  
  return(res1)
}




pv_continuous <- function(xuhao, snp, trait){
  
  snp.i <- as.numeric(snp[xuhao,])
  r1 <- table(snp.i)
  r1.na <- as.numeric(names(r1))
  r1.1 <- which(r1<8)
  if(length(r1.1)>0){
    snp.i[which(snp.i==r1.na[r1.1])] <- NA
  }
  
  geno_id <- which(!is.na(snp.i))
  pheno_id <- which(!is.na(trait))
  valid_id <- intersect(geno_id, pheno_id)
  
  geno <- snp.i[valid_id]
  pheno <- as.numeric(trait[valid_id])
  dats <- data.frame(pheno, geno)
  
  H1 <- aov(pheno ~ geno,data=dats)
  
  pvalue <- summary(H1)[[1]][[1,"Pr(>F)"]] 
  pvalue
}



net.feature <- function(gg,index){
  
  Con <- degree(gg,index,mode = "all", normalized = T)
  Connectivity <- c(mean(Con),sd(Con))
  
  Clo <- closeness(gg,index,mode = "all", normalized = T)
  Closeness <- mean(mean(Clo),sd(Clo))
  
  Bet <- betweenness(gg,index,directed =F, nobigint = TRUE, normalized = T)
  Betweenness <- c(mean(Bet),sd(Bet))
  
  Ecc <- eccentricity(gg,index, mode= "all")
  Eccentricity <- c(mean(Ecc),sd(Ecc))
  
  Eig <- evcent(gg,scale = T)$vector[index]
  Eigenvector  <- c(mean(Eig),sd(Eig)) 
  
  PR <- page.rank (gg,vids=index, damping = 0.85)$vector
  PageRank <- c(mean(PR),sd(PR)) 
  
  list(Connectivity=Con,Closeness=Clo,Betweenness=Bet,
       Eccentricity=Ecc,Eigenvector=Eig,PageRank=PR)
}




all.net <- function(res,method="ant",xuhao){
  
  ds <- res$gs
  dw <- res$gw
  
  sname <- colnames(ds)
  sname1 <- c()
  for(i in 1:length(sname)){
    sname1 <- c(sname1,strsplit(sname[i],"_S")[[1]][1])
  }
  sname2 <- unique(sname1)
  
  
  wname <- colnames(dw)
  wname1 <- c()
  for(i in 1:length(wname)){
    wname1 <- c(wname1,strsplit(wname[i],"_W")[[1]][1])
  }
  wname2 <- unique(wname1)
  
  sw <- table(c(sname2,wname2))
  swname <- names(sw)
  swname1 <- swname[which(sw==2)]
  
  ds1 <- c()
  for(i in 1:length(swname1)){
    
    index <- which(sname1==swname1[i])
    ds1 <- cbind(ds1,as.numeric(rowSums(ds[,index])))
  }
  colnames(ds1) <- swname1
  rownames(ds1) <- res$gname
  
  dw1 <- c()
  for(i in 1:length(swname1)){
    
    index <- which(wname1==swname1[i])
    dw1 <- cbind(dw1,as.numeric(rowSums(dw[,index])))
  }
  colnames(dw1) <- swname1
  rownames(dw1) <- res$gname
  
  ds1 <- ds1[,xuhao]
  dw1 <- dw1[,xuhao]
  nnn <- length(swname1)
  
  Con.s <- c();Con.w <- c();Clo.s <- c();Clo.w <- c();Bet.s <- c();Bet.w <- c();
  Ecc.s <- c();Ecc.w <- c();Eig.s <- c();Eig.w <- c();PR.s <- c();PR.w <- c();
  for(i in 1:nnn){
    s.n <- thr.compute(pheno=ds1[,i])
    w.n <- thr.compute(pheno=dw1[,i])
    s.n$ant[which(is.infinite(s.n$ant))] <- 0
    s.n$agg[which(is.infinite(s.n$agg))] <- 0
    w.n$ant[which(is.infinite(w.n$ant))] <- 0
    w.n$agg[which(is.infinite(w.n$agg))] <- 0
    
    astmp <- network.sparse(gg=s.n,en=200)
    awtmp <- network.sparse(gg=w.n,en=200)
    
    if(method=="mut"){
      g1.s <- graph.adjacency(t(astmp$mut.n),mode="directed",weighted=T)
      index1.s <- which(degree(g1.s)>10)
      g1.w <- graph.adjacency(t(awtmp$mut.n),mode="directed",weighted=T)
      index1.w <- which(degree(g1.w)>10)
    }
    
    if(method=="ant"){
      g1.s <- graph.adjacency(t(astmp$ant.n),mode="directed",weighted=T)
      index1.s <- which(degree(g1.s)>12)
      g1.w <- graph.adjacency(t(awtmp$ant.n),mode="directed",weighted=T)
      index1.w <- which(degree(g1.w)>12)
    }
    
    if(method=="agg"){
      g1.s <- graph.adjacency(t(astmp$agg.n),mode="directed",weighted=T)
      index1.s <- which(degree(g1.s)>12)
      g1.w <- graph.adjacency(t(awtmp$agg.n),mode="directed",weighted=T)
      index1.w <- which(degree(g1.w)>12)
    }
    
    if(method=="alt"){
      g1.s <- graph.adjacency(t(astmp$alt.n),mode="directed",weighted=T)
      index1.s <- which(degree(g1.s)>12)
      g1.w <- graph.adjacency(t(awtmp$alt.n),mode="directed",weighted=T)
      index1.w <- which(degree(g1.w)>12)
    }
    
    Con.s <- c(Con.s,mean(degree(g1.s,index1.s,mode = "all", normalized = T)))
    Con.w <- c(Con.w,mean(degree(g1.w,index1.w,mode = "all", normalized = T)))
    Clo.s <- c(Clo.s,mean(closeness(g1.s,index1.s,mode = "all", normalized = T)))
    Clo.w <- c(Clo.w,mean(closeness(g1.w,index1.s,mode = "all", normalized = T)))
    Bet.s <- c(Bet.s,mean(betweenness(g1.s,index1.s,directed =F, nobigint = TRUE, normalized = T)))
    Bet.w <- c(Bet.w,mean(betweenness(g1.w,index1.w,directed =F, nobigint = TRUE, normalized = T)))
    Ecc.s <- c(Ecc.s,mean(eccentricity(g1.s,index1.s,mode= "all")))
    Ecc.w <- c(Ecc.w,mean(eccentricity(g1.w,index1.w,mode= "all")))
    Eig.s <- c(Eig.s,mean(evcent(g1.s,scale = T)$vector[index1.s]))
    Eig.w <- c(Eig.w,mean(evcent(g1.w,scale = T)$vector[index1.w]))
    PR.s <- c(PR.s,mean(page.rank(g1.s, damping = 0.85)$vector[index1.s]))
    PR.w <- c(PR.w,mean(page.rank(g1.w, damping = 0.85)$vector[index1.w]))
  }
  

  return(rbind(Con.s=Con.s,Con.w=Con.w,Clo.s=Clo.s,Clo.w=Clo.w,Bet.s=Bet.s,Bet.w=Bet.w,
                 Ecc.s=Ecc.s,Ecc.w=Ecc.w,Eig.s=Eig.s,Eig.w=Eig.w,PR.s=PR.s,PR.w=PR.w))
  
}





net.cor<- function(res){
  
  allsn.mut <- c(); allwn.mut <- c();
  allsn.ant <- c(); allwn.ant <- c();
  allsn.agg <- c(); allwn.agg <- c();
  allsn.alt <- c(); allwn.alt <- c();
  for(i in 1:206){
    s.n <- thr.compute(pheno=res$gs[,i])
    w.n <- thr.compute(pheno=res$gw[,i])
    allsn.mut <- cbind(allsn.mut,s.n$mut);allwn.mut <- cbind(allwn.mut,w.n$mut)
    allsn.ant <- cbind(allsn.ant,s.n$ant);allwn.ant <- cbind(allwn.ant,w.n$ant)
    allsn.agg <- cbind(allsn.agg,s.n$agg);allwn.agg <- cbind(allwn.agg,w.n$agg)
    allsn.alt <- cbind(allsn.alt,s.n$alt);allwn.alt <- cbind(allwn.alt,w.n$alt)
  }
  nc <- dim(allsn.mut)[1]
  c.mut <- c();c.ant <- c();c.agg <- c();c.alt <- c();
  for(i in 1:nc){
    i1 <- which(is.infinite(allsn.mut[i,]))
    i2 <- which(is.infinite(allwn.mut[i,]))
    ii <- unique(c(i1,i2))
    c.mut <- c(c.mut,cor(allsn.mut[i,-ii],allwn.mut[i,-ii]))
    i1 <- which(is.infinite(allsn.ant[i,]))
    i2 <- which(is.infinite(allwn.ant[i,]))
    ii <- unique(c(i1,i2))
    c.ant <- c(c.ant,cor(allsn.ant[i,-ii],allwn.ant[i,-ii]))
    i1 <- which(is.infinite(allsn.agg[i,]))
    i2 <- which(is.infinite(allwn.agg[i,]))
    ii <- unique(c(i1,i2))
    c.agg <- c(c.agg,cor(allsn.agg[i,-ii],allwn.agg[i,-ii]))
    i1 <- which((allsn.alt[i,])==0)
    i2 <- which((allwn.alt[i,])==0)
    ii <- unique(c(i1,i2))
    c.alt <- c(c.alt,cor(allsn.alt[i,-ii],allwn.alt[i,-ii]))
  }
  mut.i1 <- which(is.infinite(c(allsn.mut)))
  mut.i2 <- which(is.infinite(c(allwn.mut)))
  mut.ii <- unique(c(mut.i1,mut.i2))
  allc.mut <- cor(c(allsn.mut)[-mut.ii],c(allwn.mut)[-mut.ii])
  
  ant.i1 <- which(is.infinite(c(allsn.ant)))
  ant.i2 <- which(is.infinite(c(allwn.ant)))
  ant.ii <- unique(c(ant.i1,ant.i2))
  allc.ant <- cor(c(allsn.ant)[-ant.ii],c(allwn.ant)[-ant.ii])
  
  agg.i1 <- which(is.infinite(c(allsn.agg)))
  agg.i2 <- which(is.infinite(c(allwn.agg)))
  agg.ii <- unique(c(agg.i1,agg.i2))
  allc.agg <- cor(c(allsn.agg)[-agg.ii],c(allwn.agg)[-agg.ii])
  
  alt.i1 <- which((c(allsn.alt))==0)
  alt.i2 <- which((c(allwn.alt))==0)
  alt.ii <- unique(c(alt.i1,alt.i2))
  allc.alt <- cor(c(allsn.alt)[-alt.ii],c(allwn.alt)[-alt.ii])
  
  list(c.mut=c.mut,c.ant=c.ant,c.agg=c.agg,c.alt=c.alt,allc.mut=allc.mut,allc.ant=allc.ant,allc.agg=allc.agg,
       allc.alt=allc.alt)
}


net.cor.bmi<- function(res,bmi){
  
  allsn.mut <- c(); allwn.mut <- c();
  allsn.ant <- c(); allwn.ant <- c();
  allsn.agg <- c(); allwn.agg <- c();
  allsn.alt <- c(); allwn.alt <- c();
  for(i in 1:206){
    s.n <- thr.compute(pheno=res$gs[,i])
    w.n <- thr.compute(pheno=res$gw[,i])
    allsn.mut <- cbind(allsn.mut,s.n$mut);allwn.mut <- cbind(allwn.mut,w.n$mut)
    allsn.ant <- cbind(allsn.ant,s.n$ant);allwn.ant <- cbind(allwn.ant,w.n$ant)
    allsn.agg <- cbind(allsn.agg,s.n$agg);allwn.agg <- cbind(allwn.agg,w.n$agg)
    allsn.alt <- cbind(allsn.alt,s.n$alt);allwn.alt <- cbind(allwn.alt,w.n$alt)
  }
  nc <- dim(allsn.mut)[1]
  cs.mut <- c();cs.ant <- c();cs.agg <- c();cs.alt <- c();
  for(i in 1:nc){
    i1 <- which(is.infinite(allsn.mut[i,]))
    cs.mut <- c(cs.mut,cor(allsn.mut[i,-i1],bmi[-i1]))
    i1 <- which(is.infinite(allsn.ant[i,]))
    cs.ant <- c(cs.ant,cor(allsn.ant[i,-i1],bmi[-i1]))
    i1 <- which(is.infinite(allsn.agg[i,]))
    cs.agg <- c(cs.agg,cor(allsn.agg[i,-i1],bmi[-i1]))
    i1 <- which((allsn.alt[i,])==0)
    cs.alt <- c(cs.alt,cor(allsn.alt[i,-i1],bmi[-i1]))
  }
  
  cw.mut <- c();cw.ant <- c();cw.agg <- c();cw.alt <- c();
  for(i in 1:nc){
    i1 <- which(is.infinite(allwn.mut[i,]))
    cw.mut <- c(cw.mut,cor(allwn.mut[i,-i1],bmi[-i1]))
    i1 <- which(is.infinite(allwn.ant[i,]))
    cw.ant <- c(cw.ant,cor(allwn.ant[i,-i1],bmi[-i1]))
    i1 <- which(is.infinite(allwn.agg[i,]))
    cw.agg <- c(cw.agg,cor(allwn.agg[i,-i1],bmi[-i1]))
    i1 <- which((allwn.alt[i,])==0)
    cw.alt <- c(cw.alt,cor(allwn.alt[i,-i1],bmi[-i1]))
  }
  list(cs.mut=cs.mut,cs.ant=cs.ant,cs.agg=cs.agg,cs.alt=cs.alt,
       cw.mut=cw.mut,cw.ant=cw.ant,cw.agg=cw.agg,cw.alt=cw.alt)
}






thr.compute <- function(pheno){
  
  n <- length(pheno)
  ncm <- combn(n,2)
  
  mut <-  apply(ncm,2,function(x,y){
    g1 <- y[x[1]]
    g2 <- y[x[2]]
    if(g1==g2){
      return(log(g1)+log(g2))
    }else
      return(log(g1)+log(g2)-log(abs(g1-g2)))
  },y=pheno)
  
  agg <-  apply(ncm,2,function(x,y){
    g1 <- y[x[1]]
    g2 <- y[x[2]]
    if(g1==g2){
      return(1)
    }else
      return(log(max(c(g1,g2)))-log(min(c(g1,g2))))
  },y=pheno)
  
  ant <-  apply(ncm,2,function(x,y){
    g1 <- y[x[1]]
    g2 <- y[x[2]]
    if(g1==g2){
      return(-(log(g1)+log(g2)))
    }else
      return(-(log(g1)+log(g2)+log(abs(g1-g2))))
    
  },y=pheno)
  
  alt <-  apply(ncm,2,function(x,y){
    g1 <- y[x[1]]
    g2 <- y[x[2]]
    if(g1==g2){
      return(0)
    }else
      return(log(max(c(g1,g2))-min(c(g1,g2)))-log(max(c(g1,g2))))
  },y=pheno)
  
  return(list(mut=mut,agg=agg,ant=ant,alt=alt,n=n,ncm=ncm))
}




network.sparse <- function(gg,en=200){
  
  ng <- gg$n
  nnc <- gg$ncm
  mut <- gg$mut
  agg <- gg$agg
  ant <- gg$ant
  alt <- gg$alt
  
  mut.n <- netwok.thre(mut,en,nnc,ng)
  agg.n <- netwok.thre(agg,en,nnc,ng)
  ant.n <- netwok.thre(abs(ant),en,nnc,ng)
  alt.n <- netwok.thre(abs(alt),en,nnc,ng)
  
  return(list(mut.n=mut.n,agg.n=agg.n,ant.n=ant.n,alt.n=alt.n))
}






netwok.thre <- function(sel,en,nnc,ng){
  
  thre <- sort(sel,decreasing = T)[en]
  mm <- matrix(0,ng,ng)
  for(ii in 1:dim(nnc)[2]){
    nnc1 <- nnc[,ii]
    mm[nnc1[1],nnc1[2]] <- as.numeric(sel[ii] > thre)
  }
  for(i in 1:(dim(mm)[1]-1)){
    i.v <- sort(sel[which(nnc[1,]==i)],decreasing = T)[1:3]
    for(j in 1:length(i.v)){
      ii.i <- nnc[,which(i.v[j]==sel)]
      mm[ii.i[1],ii.i[2]] <- 1
    }
  }
  return(mm)
}



get_con_param<-function(parm.id)
{
  for (e in commandArgs())
  {
    ta = strsplit(e,"=", fixed=TRUE);
    if(! is.na( ta[[1]][2]))
    {
      temp = ta[[1]][2];
      if( ta[[1]][1] == parm.id) {
        return (as.character(temp));
      }
    }
  }
  
  return(NA);
}


