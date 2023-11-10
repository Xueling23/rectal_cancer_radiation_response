# BiocManager::install("readxl")
# BiocManager::install("clusterProfiler")
# BiocManager::install("DOSE")

library("readxl")
library("clusterProfiler")
library("DOSE")
library(enrichplot)
library(org.Hs.eg.db)
############# load an enrichment result in one sheet #############
readD <- function(sheetNum){
  data1 <- read_xlsx("Table S1.xlsx",     # Read xlsx file with read.xlsx
                   sheet = sheetNum)
  colnames(data1)<- data1[1,]
  data1 <- data1[-1,]
  return(data1)
}

############## load all enrichment results in multiple sheets #############
loadAll <- function(){
  goBulK <- readD(3)
  goT <- readD(8)
  goB <- readD(9)
  goM <- readD(10)
  goS <- readD(11)
  goE <- readD(12)
  gA <- readD(7)
  gB <- readD(2)
  goAll <- list("goBulk" = goBulK, "goT"=goT, "goB" = goB,"goM" = goM, "goS"= goS, "goE" = goE,"gA" = gA,"gB" = gB)
  return(goAll)
}

########### Prepare enriched terms for each cell type ###############
prepD <- function(goB){
  # goB <- goA$goBulk
  # df1 <- goB[,c(1,2,4)]
  df1 <- cbind.data.frame(cbind(goB$source,goB$term_name,goB$adjusted_p_value))
  colnames(df1) <- c( "source","term_name", "adjusted_p_value")
  df1$GeneRatio <- as.numeric(goB$intersection_size)/as.numeric(goB$query_size)
  df1$Count <- as.numeric(goB$intersection_size)
  df1$fd <- df1$GeneRatio/(as.numeric(goB$term_size)/as.numeric(goB$effective_domain_size))
  
  data2 <- df1[which(df1$source=="GO:BP"),]
  # data2 <- df1[which(df1$source=="KEGG"| df1$source=="REAC"),]
  data2$P.adjust <- as.numeric(data2$adjusted_p_value)
  data2 <- data2[,-3]
  data2 <- data2[order(data2$P.adjust,decreasing = FALSE), ]
  return(data2)
}

############ enriched terms visualization ###################
datavis <- function(data2){
   #ggplot2::ggplot(data = data2[c(1:4,6:31),], ggplot2::aes(x = GeneRatio, y = term_name,
    ggplot2::ggplot(data = data2, ggplot2::aes(x = GeneRatio, y = term_name, 
                                                    color = P.adjust, size = Count)) +
    ggplot2::geom_point() +
    # ggplot2::scale_color_gradient(low ="red" , high ="blue" ) + 
    ggplot2::scale_color_gradient(low = "black", high ="grey" ) + 
    ggplot2::theme_bw() +
    ggplot2::ylab("") +
    ggplot2::xlab("") +
    ggplot2::theme(text = ggplot2::element_text(size=20) ) 
    # ggplot2::ggtitle("Enriched BP")
    # ggplot2::ggtitle("Enriched KEGG/REAC")
}
######## bulk enrichment result visualization ###########
bulkEnr <- function(goA){
  goB <- goA$goBulk
  data2 <- prepD(goB)
  # datavis(data2)
  return(data2)
}
# data2 <- prepD(goB)
# datavis(data2[c(1:4,6:31),])

######## myeloid cell type-specific enrichment visualization ###############
myeloidEnr <- function(goA){
  goB <- goA$goM
  data2 <- prepD(goB)
  # datavis(data2)
  return(data2)
}
# data2M <- myeloidEnr(goA)
# datavis(data2M)
# datavis(data2M[c(1:3,5:19,21:32),])## some extreme long terms don't show for space limitation
# datasK <- stromalEnr(goA)
# datavis(datasK)

######## stromal cell type-specific enrichment visualization ############### 
stromalEnr <- function(goA){
  goB <- goA$goS
  data2 <- prepD(goB)
   return(data2)
}
# data2S <- stromalEnr(goA)
# datavis(data2S[c(1:12,14:31),]) ## some extreme long terms don't show for space limitation

# saveRDS(data2B,file = "GO_bulk.rds")
# saveRDS(data2S,file = "GO_stromal.rds")
# saveRDS(data2M,file = "GO_myeloid.rds")
# saveRDS(dataBK,file = "KEGG_RECT_bulk.rds")
# saveRDS(datasK,file = "KEGG_RECT_stromal.rds")
# saveRDS(dataMK,file = "KEGG_RECT_myeloid.rds")

############# DOSE and clusterProfiler based visualization methods ##############
######## subtle difference may exist due to different enrichment analysis tools
goVis <- function(gB){
  de <- gB$ID
  de <- bitr(de,fromType ='SYMBOL' ,toType = 'ENTREZID',OrgDb = "org.Hs.eg.db")
  ego <- enrichGO(de$ENTREZID, OrgDb = "org.Hs.eg.db", ont="BP", readable=TRUE)
  dotplot(ego,showCategory=20,font.size = 14)
  return(de)
}

keggVis <- function(de){
  kk <- enrichKEGG(de$ENTREZID, organism = "hsa",  pvalueCutoff = 0.05,qvalueCutoff = 0.05)
  p <- dotplot(kk,showCategory=20)
}
