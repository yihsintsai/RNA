library(edgeR)
library("RUVSeq")
library(data.table)
library(RColorBrewer)
library(dplyr)
library("ggcorrplot")
source("/mnt/nas2/yh/10.RNA_seq_poly_A/3.R/ercc_cor.R")
setwd("/mnt/nas2/yh/10.RNA_seq_poly_A/8.count_ucsc/")

sample_list <- data.frame(read.table("sample_list.txt",col.names ="sampleId" ))
gene_list <-read.table("gene_id.list",col.names=c("geneID"))
   extra_list<-c(10,11,12,16,28,32)
data <- data.frame(gene_list)
#for (f in 1:36){
for (f in 1:nrow(sample_list)){
  sample <- sample_list[f,]
  sample_txt <- paste0(sample,"_R1.fq.gz.HQ_filter_count_for_r.txt",sep = "")
  table <- data.frame(read.table(sample_txt,col.names=c("geneID",sample)))%>%
    dplyr::select(2)
  data <- cbind(data,table) 
}

row.names(data) <- data$geneID

RNA.table <- data[,2:51]
    
    rna_wgs_replace_name <-  data.frame(read.table("/mnt/nas2/yh/10.RNA_seq_poly_A/SNP/rna_wga_samplename_op", header=FALSE))
    rna.replace <- rna_wgs_replace_name[grepl("rna",rna_wgs_replace_name$V1),]
    rna.replace <- sub("rna_", "", rna.replace)

row.count <- data.frame(RNA.table)
    colnames(row.count)<-rna.replace
    row.count <- row.count[, order(colnames(row.count))]
    row.count <- row.count[,c(30:50,1:29)]


####Filter data #######################################################################
#requiring > 5 reads in at least two samples for each gene
filter <- apply(row.count , 1, function(x) length(x[x>5])>=2)
    filtered <- row.count[filter,]

##ERCC data prepare
ERCC <- rownames(filtered)[grep("^ERCC-", rownames(filtered))]
    ERCC_count <- filtered[row.names(filtered)%like%"^ERCC-",]


#ERCC_correlation 
ERCC_con_table <- read.csv("/mnt/nas2/yh/10.RNA_seq_poly_A/reference/ERCC_ref/ERCC_control_table.csv") %>% 
  dplyr::select(1:4) %>% 
  `names<-`(c("RE-sort_ID","ERCC_ID","subgroup","Conc.Mix1"))


condition <- as.factor(rep(c("SHORT","LONG","SHORT","LONG","SHORT","LONG"),c(19,17,1,6,1,6)))  #short/long
    ERCC_correlation <- ercc_cor(ercc_control.table = ERCC_con_table,rna_table = RNA.table,condition = condition,corr_arg = "spearman")

##RNA sample data prepare 
genes <- filtered[!grepl("^ERCC-",row.names(filtered)),]

##datacondition set
data <- newSeqExpressionSet(as.matrix(filtered),
                            phenoData = data.frame(condition, row.names=colnames(filtered)))



####plot non_nomorlization#######################################################################
##plot non_nomorlization
colors <- brewer.pal(3, "Dark2")
plotRLE(data, outline=FALSE, ylim=c(-4, 4), col=colors[condition])
plotPCA(data, col=colors[condition], cex=0.8)

##plot bwt_nomorlization
bwt_data <- betweenLaneNormalization(data, which="upper")
    plotRLE(bwt_data, outline=FALSE, ylim=c(-4, 4), col=colors[condition])
    plotPCA(bwt_data, col=colors[condition], cex=0.8)

##plot ERCC_nomorlization
ercc_data <- RUVg(bwt_data, ERCC, k=3)
    plotRLE(ercc_data, outline=F, ylim=c(-4, 4), col=colors[condition],
        main = 'RNA-seq without unwant variation k=5')
    plotPCA(ercc_data, col=colors[condition], cex=0.8,
        main = 'RNA-seq without unwant variation k=5')

rna_count.data <- data.frame(ercc_data@assayData[["normalizedCounts"]])
    summary(rna_count.data)



####count cpm#######################################################################
##count cpm
#normalize with ERCC
gene_length <- read.table("/mnt/nas2/yh/10.RNA_seq_poly_A/reference/htseq/ucsc/bed/chm13.draft_v2.0.gene_annotation.sorted.filtered.vfilter_ky_ercc.gene_length", header=F) %>%
  `colnames<-`(c("name","gene.length"))
gene_name<- data.frame(rownames(rna_count.data))%>%`colnames<-`(c("name"))
gene_length.table <- merge(gene_name, gene_length, by="name")%>%
  `colnames<-`(c("GeneID","Lengths"))
#gene_length_info <- data.frame(GeneID = gene_ids, Length = gene_lengths)
design <- model.matrix(~condition + W_1 + W_2 + W_3 ,  data = pData(ercc_data))
    y <- DGEList(counts=counts(ercc_data), group=condition ,genes=gene_length.table)
    y <- calcNormFactors(y, method="upperquartile")
    y <- estimateGLMCommonDisp(y, design)
    y <- estimateGLMTagwiseDisp(y, design)

    normalize_cpm.table <- data.frame(cpm(y))
    normalize_rpkm.table <- data.frame(rpkm(y))


####filter XY_gene#######################################################################

XYM_gene <- read.table("/mnt/nas2/yh/10.RNA_seq_poly_A/reference/htseq/ucsc/XYM_gene.list")

normalize_cpm.table <- normalize_cpm.table[!(rownames(normalize_cpm.table) %in% XYM_gene$V1), ]


########################################################
## merge with oncogene
OncogeneList <- read.csv("../OncogeneList.csv")
final_onco.table <- final_count.table[row.names(final_count.table)%in%OncogeneList$Gene.Symbol,]
write.csv(final_onco.table,"../final_onco.table.csv")
result_onco_data <- result_data[row.names(result_data)%in%OncogeneList$Gene.Symbol,]
write.csv(result_onco_data ,"../fresult_onco_data.csv")
gain_lose_count <- final_count.table[row.names(final_count.table)%in%gain_lose$symbol,]
write.csv(gain_lose_count ,"../gain_lose_count.csv")
gain_lose_FC <- result_data[row.names(result_data)%in%gain_lose$symbol,]
write.csv(gain_lose_FC ,"../gain_lose_FC.csv")

########################################################
## plot spike normalized count heatmap
 #corrplot
#pearson
pearson <-cor(final_count.table,method = 'pearson')
round(pearson, 2)
col<- colorRampPalette(c("#7744FF", "white","#F50A0A"))(299)
corrplot(pearson, 
         method="color",
         col=col, 
         number.cex = 0.8, 
         col.lim = c(min(0.75), max(1)),
         tl.cex=1.1,
         type="full", 
         order="hclust", 
         addCoef.col = "black",  # Add coefficient of correlation
         tl.col="black",
         tl.srt=90,  #Text label color and rotation
         is.corr = FALSE)
#spearman
spearman <- cor(final_count.table,method = 'spearman')
round(spearman, 2)
corrplot(spearman, 
         method="color",
         col=col, 
         number.cex = 0.8, 
         col.lim = c(min(0.75), max(1)),
         tl.cex=1.1,
         type="full", 
         order="hclust", 
         addCoef.col = "black",  # Add coefficient of correlation
         tl.col="black",
         tl.srt=90,  #Text label color and rotation
         is.corr = FALSE)
#########################################################
#SHORT/LONG heatmap
short <- data.frame(final_count.table[,1:19])
short_res3 <- cor(short,method = 'pearson')
round(short_res3, 2)
col_2<- colorRampPalette(c("white","#F58484","#F50A0A"))(299)
corrplot(short_res3, 
         method="color", 
         col=col_2, 
         number.cex = 1.3,
         tl.cex=1.5 ,
         col.lim = c(min(0.75), max(1)),
         type="full", 
         order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=90,#Text label color and rotation
         is.corr = FALSE)


long <- data.frame(final_count.table[,20:36])
long_res3 <- cor(long,method = 'pearson')
round(long_res3, 2)
col_3<- colorRampPalette(c("white","#CAB1FA","#7744FF"))(299)
corrplot(long_res3, method="color", col=col_3, number.cex = 1.3, tl.cex=1.5,
         col.lim = c(min(0.75), max(1)),
         type="full", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=90,#Text label color and rotation
         is.corr = FALSE)

#################################################
##copmare condition
#vocano plot
#vocanal filter
up_gene <- subset(result_data, result_data$logFC > 0 & result_data$PValue < 0.05)
#up : 385
up5_gene <-subset(up_gene, up_gene$logFC > 4.7 | up_gene$PValue < 0.0000005)

down_gene <- subset(result_data, result_data$logFC < -1 & result_data$PValue < 0.05)
#down : 147
down5_gene <- subset(down_gene, down_gene$logFC < -1 | down_gene$PValue < 0.0000005)

#merge oncogene
up.onco.up.gene <- up_gene[row.names(up_gene)%in%OncogeneList$Gene.Symbol,]
write.csv(up.onco.up.gene,"../up_oncogene.csv")
down.onco.down.gene <- down_gene[row.names(down_gene)%in%gain_lose$symbol,]
write.csv(down.onco.down.gene,"../down_oncogene.csv")

### plot
plot(
  result_data$logFC, -log2(result_data$PValue),
  pch = 19, col = "grey", cex = 0.5, xlim = c(-10,10), ylim = c(0,20),
  ylab = " ", xlab = " ", xaxt = "n", yaxt = "n")
par(new=TRUE)
plot(
  up_gene$logFC, -log2(up_gene$PValue),
  pch = 19, col = "red", cex = 0.5, xlim = c(-10,10), ylim = c(0,20),
  ylab = " ", xlab = " ", xaxt = "n", yaxt = "n")
par(new=TRUE)
plot(
  down_gene$logFC, -log2(down_gene$PValue),
  pch = 19, col = "blue", cex = 0.5, xlim = c(-10,10), ylim = c(0,20),
  ylab = "-Log2(Pval)", xlab = "Log2FC (short vs. long)")
abline(v=1); abline(v=-1); abline(h=-log2(0.05))
par(new=TRUE)
plot(
  down.onco.down.gene$logFC, -log2(down.onco.down.gene$PValue),
  pch = 19, col = "green", cex = 0.5, xlim = c(-10,10), ylim = c(0,20),
  ylab = "-Log2(Pval)", xlab = "Log2FC (short vs. long)")
text(down.onco.down.gene$logFC, -log2(down.onco.down.gene$PValue),#row.names(down.onco.down.gene),
     cex=1, pos=3,col="black") 
#text(up5_gene$logFC, -log2(up5_gene$PValue),row.names(up5_gene),
     #cex=0.7, pos=3,col="black") 
#text(down5_gene$logFC, -log2(down5_gene$PValue),row.names(down5_gene),
     #cex=0.7, pos=3,col="black") 
#text(up_gene_onco$logFC, -log2(up_gene_onco$PValue),row.names(up_gene_onco),
     #cex=2, pos=3,col="black") 
#text(down_gene_onco$logFC, -log2(down_gene_onco$PValue),row.names(down_gene_onco),
     #cex=0.7, pos=3,col="black")

###########################################################
###gene ontology
#gene GO
library(clusterProfiler); library(org.Hs.eg.db); library(dplyr) 
library(enrichplot)
up.GO <- row.names(up_gene) %>% enrichGO(
  gene=., OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "ALL")

up_GO_data <- data.frame(up.GO)
write.csv(up_GO_data ,"../up_GO_data.csv")

dotplot(up.GO, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free") #+
  #scale_colour_gradient2(low = "#FF00FF", mid = "#000000", high = "#FFFF00")

dotplot(up.GO, showCategory = 13)
up_GO_BP <- up.GO[up.GO@result[["ONTOLOGY"]] %in% "BP",]
up_GO_BP <- up_GO_BP[order(up_GO_BP$p.adjust,decreasing = T),]

library(ggplot2) 
library(dplyr) 
showCategory =15
font.size =12
p<-up_GO_BP %>% 
  slice(1:showCategory) %>% 
  ggplot(aes(x=forcats::fct_reorder(Description,p.adjust,.desc = T),y=Count,fill=p.adjust))+ 
  geom_bar(stat="identity")+
  coord_flip()+
  scale_fill_continuous(low="red", high="blue", guide=guide_colorbar(reverse=TRUE))+
  labs(x=NULL) +
  ggtitle("")+
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black",
                                   size = font.size, vjust =1 ),
        axis.text.y = element_text(colour = "black",
                                   size = font.size, hjust =1 ),
        axis.title = element_text(margin=margin(10, 5, 0, 0),
                                  color = "black",size = font.size),
        axis.title.y = element_text(angle=90))
p
#################################################################################
#up.onco.up.gene
up.onco.up.gene <- row.names(up.onco.up.gene) %>% enrichGO(
  gene=., OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "ALL")

up.onco.up.gene_BP  <- up.onco.up.gene [up.onco.up.gene@result[["ONTOLOGY"]] %in% "BP",]
up.onco.up.gene_BP  <- up.onco.up.gene_BP[order(up.onco.up.gene_BP $p.adjust,decreasing = F),]
up.onco.up.gene_BP <- up.onco.up.gene_BP[,1:10]
library(ggplot2) 
library(dplyr)
showCategory =15
font.size =12
p<-up.onco.up.gene_BP %>% 
  slice(1:showCategory) %>% 
  ggplot(aes(x=forcats::fct_reorder(Description,p.adjust),y=Count,fill=p.adjust))+ 
  geom_bar(stat="identity")+
  coord_flip()+
  scale_fill_continuous(low="red", high="blue", guide=guide_colorbar(reverse=TRUE))+
  labs(x=NULL) +
  ggtitle("")+
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black",
                                   size = font.size, vjust =1 ),
        axis.text.y = element_text(colour = "black",
                                   size = font.size, hjust =1 ),
        axis.title = element_text(margin=margin(10, 5, 0, 0),
                                  color = "black",size = font.size),
        axis.title.y = element_text(angle=90))
p



