# hnsc <- read.table("D:/工作/导师/建工/文章/LXL/Rejected/重修/返回/TCGA-HNSC.htseq_counts.tsv/TCGA-HNSC.htseq_counts.tsv",sep = "\t",header=T,check.names=F,row.names=1)
luad_ph <- read.table("D:/工作/导师/建工/文章/LXL/Rejected/重修/返回/LUAD_TCGA/Phenotype_LUAD.TXT",sep = "\t",header=T,check.names=F,row.names=1)
luad_exp <- read.table("D:/工作/导师/建工/文章/LXL/Rejected/重修/返回/TCGA-LUAD.htseq_counts.tsv/TCGA-LUAD.htseq_counts.tsv",sep = "\t",header=T,check.names=F,row.names=1)
luad_ph_r <- read.csv("I:/培养/TCGA-LUAD.GDC_phenotype.csv",sep = ",",header=T,check.names=F,row.names=1)
luad_ph1 <- luad_ph_r[c("additional_radiation_therapy","followup_treatment_success")]
luad_ph2 <- luad_ph1[luad_ph1$additional_radiation_therapy=="YES",]
luad_ph3 <- luad_ph2[luad_ph2$followup_treatment_success!="",]
luad_ph4 <- luad_ph3[intersect(colnames(luad_exp),row.names(luad_ph3)),]
luad_exp2 <- luad_exp[intersect(colnames(luad_exp),row.names(luad_ph3))]
library(org.Hs.eg.db)
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = luad_exp2,
                              colData = luad_ph4,
                              design = ~ nr)
dds <- DESeq(dds)
res <- results(dds)
sigG <- row.names(res)[which(res$pvalue<0.01)]
sigG_r <- str_split(sigG,"\\.")
sigG_r2 <- vector()
for(i in 1:length(sigG_r)){sigG_r2[i] <- sigG_r[[i]][1]}
write.csv(de,"LUAD_bulk_DEGseq2.csv")
de <- bitr(sigG_r2,fromType = 'ENSEMBL' ,toType = c('SYMBOL','ENTREZID'),OrgDb = "org.Hs.eg.db")
de <- bitr(de$ENTREZID,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = "org.Hs.eg.db")
#csgene_EN <- bitr(csDEGs,fromType = 'SYMBOL',toType = 'ENSEMBLE',OrgDb = "org.Hs.eg.db")
go <- enrichGO(de$ENTREZID, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
close