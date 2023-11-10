grpTest <- function(geps1,geps2,stderr1,stderr2,thres){
   setwd("D:\\工作\\导师\\建工\\文章\\LXL\\20220831组模式差异分析(项目打包)\\20220831组模式的差异分析\\datasets\\GSE119409")
  
  # geps1=geps$`119409_R`
  # geps2=geps$`119409_NR`
  # stderr1 = geps$`119409_R_std`
  # stderr2 = geps$`119409_NR_std`
  qval <-list()
  zval <- list()
  q_z <-list()
  gnames <- list()
  mean_R_f <- list()
  mean_NR_f <- list()

 
   for (i in 1:ncol(geps1)){

    #### r1 and r2 represent average expressions of the gene j in cell type i in two groups respectively   
    r1 <- geps1[,i][!is.na(geps1[,i]) & geps1[,i]!=1 & !is.na(geps2[,i]) & geps2[,i]!=1]
    g1 <- row.names(geps1)[!is.na(geps1[,i]) & geps1[,i]!=1 & !is.na(geps2[,i]) & geps2[,i]!=1] 
    r2 <- geps2[,i][!is.na(geps1[,i]) & geps1[,i]!=1 & !is.na(geps2[,i]) & geps2[,i]!=1]
    g2 <- row.names(geps2)[!is.na(geps1[,i]) & geps1[,i]!=1 & !is.na(geps2[,i]) & geps2[,i]!=1] 
    r1 =as.data.frame(r1)
    r2 = as.data.frame(r2) 
    row.names(r1) <- g1
    row.names(r2) <- g2
    r <- merge.data.frame(r1,r2,by=0)
    row.names(r)<-r[,1]
    r<-r[,-1]
    colnames(r)<-c("R","NR")
    sd1 <- stderr1[,i][!is.na(geps1[,i]) & geps1[,i]!=1 & !is.na(geps2[,i]) & geps2[,i]!=1]
    sd2 <- stderr2[,i][!is.na(geps1[,i]) & geps1[,i]!=1 & !is.na(geps2[,i]) & geps2[,i]!=1]
    vBetaZ <- sapply(1:nrow(r), function(j) (r$R[j]-r$NR[j])/sqrt(sd1[j]^2+sd2[j]^2))
    ZPS <- 2*pnorm(-abs(vBetaZ))
    Zqvals <- p.adjust(ZPS, method = "BH")
    
    zval[[i]] =vBetaZ[which(Zqvals < thres)]
    gnames[[i]] =row.names(r)[which( Zqvals < thres)]
    mean_R_f[[i]] = r$R[which(Zqvals < thres)]
    mean_NR_f[[i]] = r$NR[which(Zqvals < thres)]
    qval[[i]] = Zqvals[which(Zqvals < thres)]
    
    gnames[[i]] = gnames[[i]][is.finite(qz3493_005[[i]]$Zval)]
    mean_R_f[[i]] =mean_R_f[[i]][is.finite(qz3493_005[[i]]$Zval)]
    mean_NR_f[[i]] = mean_NR_f[[i]][is.finite(qz3493_005[[i]]$Zval)]
    qval[[i]] =qval[[i]][is.finite(qz3493_005[[i]]$Zval)]
    zval[[i]] =zval[[i]][is.finite(qz3493_005[[i]]$Zval)]
      
    q_z[[i]]=list(Genes = gnames[[i]],Qval = qval[[i]],Zval =zval[[i]],Mean_R=mean_R_f[[i]],Mean_NR=mean_NR_f[[i]])
    
  }
  names(q_z) <- colnames(geps1)
  return(q_z)
  
}