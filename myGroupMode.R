#################### group-mode deconvolution ##################### 
CIBERSORTg <- function(bulk_file_txt,signature_file_txt){

 setwd("D:\\工作\\导师\\安医\\FXY\\新建文件夹\\投稿\\大修\\验证集")
 bulk_file <- read.table("Bulk expression.TXT",header=T,sep="\t",check.names=F,row.names=1)
 sigature_file <- read.table("SIG_GSE169147_2.TXT",header=T,sep="\t",check.names=F,row.names=1)
 source("D:/工作/导师/安医/FXY/新建文件夹/CIBERSORT.R")
  
  ### Note: Newman's CIBERSORTx is patented and has an extra bootstrapping procedure, variance filtering, which were not provided here,
  # although all of which were adopted in the main results of our paper. 
  
  # 
  # setwd("Path to bulk expression and signature fies")
  # bulk_file <- read.table("Bulk expression.TXT",header=T,sep="\t",check.names=F,row.names=1)
  # sigature_file <- read.table("Cell_type_signature.TXT",header=T,sep="\t",check.names=F,row.names=1)
  # source("CIBERSORT.R")
  

  c_group <- dataAugm(tpm[,c(1,4,5)])
  c_result <- CIBERSORT(X, c_group, perm=0, QN=FALSE)
  bn_group <- dataAugm(tpm[,c(2,3,6:9)]) # only cell types not equal to zeros are considered.
  bn_result <- CIBERSORT(X, bn_group, perm=0, QN=FALSE)
  
  bn_fit <- nnlsfit(bn_result,bn_group)
  c_fit <- nnlsfit(c_result,c_group)
  
  vBetaZ <- calZ (bn_fit,c_fit)
  g_sig <- gtSig(bn_group,vBetaZ)
  return(g_sig)
}
#########  nnls fitting ####################
nnlsfit <-function(c_result,c_group){
  
  # c_result <- ctr_expr_result$c_result
  # c_group <- ctr_expr_result$c_group
  # c_result <- c_result[,c(2:6,8)]
  c_result <- c_result[,c(1:2,4:6,8)] ## cell types of not-all-zero cell fractions of burn dataset.
  nr <- dim(c_group)[1] # number of genes
  x_fit <- data.frame(matrix(data=NA, nrow = nr ,ncol=6))
  p_all <- data.frame(matrix(data=NA, nrow = nr ,ncol=6))
  se_fit <- data.frame(matrix(data=NA, nrow = nr ,ncol=6))
  for(i in 1:dim(c_group)[1]){  
    mod1 <- nnls(as.matrix(c_result),as.numeric(c_group[i,]))
    x_fit[i,] <- mod1$x  #  mod1$x = coefficients(mod1) = coef(mod1)
    se_fit[i,] <-cal_t(mod1)   ###### calulate T values by executing cal_t function below.
    t1 <- as.numeric(mod1$x/se_fit[i,])
    df1 <- dim(c_group)[2]-dim(c_result)[2]-1
    # p_val <- 2*pt(q=t1, df=63, lower.tail=FALSE) # df = 70-6-1 =63
    p_val <- 2*pt(q=t1, df=df1, lower.tail=FALSE)
    p_all[i,] <- p_val
  }

  fdrs <- data.frame(matrix(data=NA, nrow = nr ,ncol=6))
  for(i in 1:6){
  fdrs[,i]<-p.adjust (p_all[,i], method="BH")
  }
  result <- list(x_fit = x_fit, p_all = p_all, fdr = fdrs,se_fit = se_fit)
}
####### T value calculation with cal_t function ########## 
cal_t<- function(mod1){
  SS_E <- t(mod1$residuals)%*%(mod1$residuals)    #SS_E = t(Y-Y_hat)%*%(Y-Y_hat)
  MS_E <- SS_E/63                   #df_E = n-p-1 = 70-6-1 = 63
  X <- as.matrix(c_result)
  # Y <- as.numeric(c_group[1,])
  inv_XtX=solve(t(X)%*%X)
  C_ii <- rep(0,6)
  se_t <- rep(0,6)
  x <- mod1$x
  for(i in 1:6){
    C_ii [i]<- inv_XtX[i,i]
    se_t[i] <- sqrt(C_ii[i]*MS_E)
   # t1[i] <- x[i]/sqrt(C_ii[i]*MS_E)
  }
  return(se_t)
}
############## Z calculation ##################
calZ <- function(bn_fit,c_fit){
  vBetaZ <- list()
  vBetaZ[[1]] <- (bn_fit$x_fit[,2]-c_fit$x_fit[,1])/sqrt(bn_fit$se_fit[,2]^2+c_fit$se_fit[,1]^2)
  for (i in 2:5){
    vBetaZ[[i]] <- (bn_fit$x_fit[,i+1]-c_fit$x_fit[,i+1])/sqrt(bn_fit$se_fit[,i+1]^2+c_fit$se_fit[,i+1]^2)
  }
  return(vBetaZ)
  } 
  

######## extract significant genes ###############

gtSig <- function(bn_group,vBetaZ){
  qvals <- list()
  vBetaz_sig <-list()
  # vBetaZ_a <- list()
  ix_sig2 <- list()
  gene_sig <- list()
  allg <- row.names(bn_group)
  for (i in 1:5){
  # vBetaZ <- (bn_fit$x_fit[,i]-c_fit$x_fit[,i])/sqrt(bn_fit$se_fit^2+c_fit$se_fit^2)
    ZPS <- 2*pnorm(-abs(as.matrix(vBetaZ[[i]])))
    qvals[[i]] <- p.adjust(ZPS, method = "BH")
    ix_sig2[[i]] <- which(bn_fit$p_all[,i]<0.05 & c_fit$p_all[,i]<0.05 & qvals[[i]]<0.05)  ## option 1, this is stringent option.
    # ix_sig2[[i]] <- which(qvals[[i]]<0.05)  ## option 2 is less stringent.
    gene_sig [[i]] <- allg[ix_sig2[[i]]]
    vBetaz_sig[[i]] <- vBetaZ[[i]][[1]][ix_sig2[[i]]]
    # vBetaZ_a[[i]] <- vBetaZ
  }
  names(ix_sig2) <- c(  "Fibroblasts","Macrophage","Keratinocytes", "CD4_T",  "Endothelial")
  names(gene_sig) <- c(  "Fibroblasts","Macrophage","Keratinocytes", "CD4_T", "Endothelial")
  names(vBetaz_sig) <- c(  "Fibroblasts","Macrophage","Keratinocytes", "CD4_T", "Endothelial")

  g_sig <- list(ix_sig = ix_sig2,gene_sig = gene_sig, vBetaz_sig = vBetaz_sig )
  return(g_sig)
}
############
# ix_sig <- list()
# ix_sig[[1]] <- which(bn_fit$p_all[,2]<0.05&c_fit$p_all[,1]<0.05)
# for(i in 2:5){
#   ix_sig[[i]] <- which(bn_fit$p_all[,i+1]<0.05&c_fit$p_all[,i+1]<0.05)
# }
