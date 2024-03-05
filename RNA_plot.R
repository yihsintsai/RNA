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

rna_wgs_replace_name <-  data.frame(read.table("/mnt/nas2/yh/10.RNA_seq_poly_A/SNP/rna_wga_samplename_op",header=FALSE))
  rna.replace <- rna_wgs_replace_name[grepl("rna",rna_wgs_replace_name$V1),]
  rna.replace <- sub("rna_", "", rna.replace)

row.count <- data.frame(RNA.table)
  colnames(row.count)<-rna.replace
  row.count <- row.count[, order(colnames(row.count))]
  row.count <- row.count[,c(30:50,1:29)]

  #L32 <-   row.count[,grepl("L32",colnames(row.count))]
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


condition <- as.factor(rep(c("SHORT","LONG"),c(21,29)))  #short/long
  ERCC_correlation <- ercc_cor(ercc_control.table = ERCC_con_table,rna_table = row.count,condition = condition,corr_arg = "spearman")
##RNA sample data prepare 

genes <- filtered[!grepl("^ERCC-",row.names(filtered)),]

##datacondition set
data <- newSeqExpressionSet(as.matrix(filtered),
                            phenoData = data.frame(condition, row.names=colnames(filtered)))



####plot non_nomorlization#######################################################################
##plot non_nomorlization
colors <- brewer.pal(3, "Dark2")
  plotRLE(data, outline=FALSE, ylim=c(-4, 4), col=colors[condition],las=2,
          main = 'RNA-seq un-normalized data')
  plotPCA(data, col=colors[condition], cex=0.8,
          main = 'RNA-seq un-normalized data')

##plot bwt_nomorlization
bwt_data <- betweenLaneNormalization(data, which="upper")
  plotRLE(bwt_data, outline=FALSE, ylim=c(-4, 4), col=colors[condition],las=2)
  plotPCA(bwt_data, col=colors[condition], cex=0.8)

##plot ERCC_nomorlization
 # for (f in 1:5) {
 #   print(f)
 #   ercc_data <- RUVg(bwt_data, ERCC, k=1)
 #     num_name = paste('RNA-seq without unwant variation k=',f, collapse = '')
 #   plotRLE(ercc_data, outline=F, ylim=c(-4, 4), col=colors[condition],
 #           main = num_name)
 #   plotPCA(ercc_data, col=colors[condition], cex=0.8,
 #           main = num_name)
 #   }
  ercc_data <- RUVg(bwt_data, ERCC, k=3)
    plotRLE(ercc_data, outline=F, ylim=c(-4, 4), col=colors[condition],
        main = 'RNA-seq without unwant variation k=3',las =3)
plotPCA(ercc_data, col=colors[condition], cex=0.8,
        main = 'RNA-seq without unwant variation k=3')

rna_count.data <- data.frame(ercc_data@assayData[["normalizedCounts"]])
  summary(rna_count.data)



####count cpm#######################################################################
##count cpm
#normalize with ERCC
gene_length <- read.table("/mnt/nas2/yh/10.RNA_seq_poly_A/reference/htseq/ucsc/bed/chm13.draft_v2.0.gene_annotation.sorted.filtered.vfilter_ky_ercc.gene_length", header=F) %>%
                  `colnames<-`(c("name","gene.length"))
  
gene_name<- data.frame(rownames(rna_count.data))%>%
                  `colnames<-`(c("name"))
gene_length.table <- merge(gene_name, gene_length, by="name")%>%
                  `colnames<-`(c("GeneID","Lengths"))
#gene_length_info <- data.frame(GeneID = gene_ids, Length = gene_lengths)
design <- model.matrix(~ W_1+W_2+W_3  + condition,  data = pData(ercc_data))
y <- DGEList(counts=counts(ercc_data), group=condition ,genes=gene_length.table)
  y <- calcNormFactors(y, method="upperquartile")
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)

normalize_cpm.table <- data.frame(cpm(y))
normalize_rpkm.table <- data.frame(rpkm(y))


####filter XY_gene#######################################################################

XYM_gene <- read.table("/mnt/nas2/yh/10.RNA_seq_poly_A/reference/htseq/ucsc/XYM_gene.list")

normalize_cpm.table <- normalize_cpm.table[!(rownames(normalize_cpm.table) %in% XYM_gene$V1), ]
 

####change sample name####################################################################
rna_wgs_replace_name <-  data.frame(read.table("/mnt/nas2/yh/10.RNA_seq_poly_A/SNP/rna_wga_samplename_op",header=FALSE))

rna.replace <- rna_wgs_replace_name[grepl("rna",rna_wgs_replace_name$V1),]
rna.replace <- sub("rna_", "", rna.replace)

row.count_filter <- data.frame(row.count)
  row.count_filter <- row.count_filter[!(rownames(row.count_filter) %in% XYM_gene$V1), ]
  write.table(row.count_filter,"/mnt/nas2/yh/20240105_rna_rowcount.table")

cpm <-  data.frame(normalize_cpm.table)

    
  
normalize_cpm.table_long <-  cpm[,grep("L", colnames(cpm))]
  normalize_cpm.table_long <- normalize_cpm.table_long[, order(colnames(normalize_cpm.table_long))]
   
normalize_cpm.table_short <-  cpm[,grep("S", colnames(cpm))]
  normalize_cpm.table_short <- normalize_cpm.table_short[, order(colnames(normalize_cpm.table_short))]
  
corre.table_short <- cor(normalize_cpm.table_short,method = 'spearman')
corre.table_long <- cor(normalize_cpm.table_long,method = 'spearman')
corre.table <- cor(cpm ,method = 'spearman')

round(corre.table_short, 2)
round(corre.table_long, 2)
round(corre.table, 2)


#####
fig <-pheatmap(corre.table_short ,
                   #clustering_method = "ward.D2",
                   color = colorRampPalette(c("white", "#FF5151"))(299),
                   cex=1.5,
                   border_color = FALSE,
                   #clustering_distance_cols="euclidean",
                   #annotation_col = annotation_col,
                   #annotation_row = annotation_col,
                   #annotation_colors = col,
                   show_rownames=T,
                   show_colnames=T,
                   cluster_row =F,
                   cluster_col =F ,
                   cex.main=0.1,
                   display_numbers = TRUE,
)
               
fig <-pheatmap(corre.table_long ,
                   #clustering_method = "ward.D2",
                   color = colorRampPalette(c("white", "#FF5151"))(299),
                   cex=1,
                   border_color = FALSE,
                   #clustering_distance_cols="euclidean",
                   #annotation_col = annotation_col,
                   #annotation_row = annotation_col,
                   #annotation_colors = col,
                   show_rownames=T,
                   show_colnames=T,
                   cluster_row =F,
                   cluster_col =F ,
                   cex.main=0.1,
                   display_numbers = TRUE
)

fig <-pheatmap(corre.table ,
               #clustering_method = "ward.D2",
               color = colorRampPalette(c("white","#FF7575", "#FF5151"))(299),
               cex=1,
               border_color = FALSE,
               #clustering_distance_cols="euclidean",
               #annotation_col = annotation_col,
               #annotation_row = annotation_col,
               #annotation_colors = col,
               show_rownames=T,
               show_colnames=T,
               cluster_row =T,
               cluster_col =T ,
               cex.main=0.1,
               display_numbers = TRUE
)
######
short <- ggcorrplot(corre.table_short,lab = TRUE) +
  scale_fill_gradient2(low = "blue", high = "red",midpoint = 0.5, breaks=c(0.5, 1), limit=c(0.5, 1))+
  theme(
    text = element_text(size = 12 , color = "black"),  
    axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1), 
    legend.position = c(1.03,0.94),
    legend.title = element_blank()
  )
plot(short)


long <- ggcorrplot(corre.table_long, lab = TRUE) +
  scale_fill_gradient2(low = "blue", high = "red", midpoint = 0.76, breaks=c(0.75, 1), limit=c(0.75, 1))+
  theme(
    text = element_text(size = 12 , color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
    legend.position = c(1.03,0.94),
    legend.title = element_blank()
  )
plot(long)


all <- ggcorrplot(corre.table,lab = TRUE) +
  scale_fill_gradient2(low = "white", high = "red",midpoint = 0.76, breaks=c(0.75, 1), limit=c(0.75, 1))+
  theme(
    text = element_text(size = 20 , color = "black"), 
    axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1), 
    legend.position = c(1.03,0.97),
    legend.title = element_blank()
  )
plot(all)
#########################################
#APOBEC list
APOBEC <- read.table("/mnt/nas2/yh/10.RNA_seq_poly_A/reference/APOBEC.list")%>%`colnames<-`(c("gene_id","gene_name"))
APOBECA.table <- data.frame(t(normalize_cpm.table[rownames(normalize_cpm.table)%in%"CHM13_G0037458",]))
  APOBECA.table$sample <- rownames(APOBECA.table)
  APOBECA.table <- APOBECA.table[order(APOBECA.table$CHM13_G0037458,decreasing =F),]
APOBECB.table <- data.frame(t(normalize_cpm.table[rownames(normalize_cpm.table)%in%"CHM13_G0037459",]))
  APOBECB.table$sample <- rownames(APOBECB.table)
  APOBECB.table <- APOBECB.table[order(APOBECB.table$CHM13_G0037459,decreasing =F),]

APOBEC.table <- normalize_cpm.table[rownames(normalize_cpm.table)%in%APOBEC$gene_id,]
APOBEC.table$gene_id <- rownames(APOBEC.table)
APOBEC.table <- merge(APOBEC.table, APOBEC, by = "gene_id")
rownames(APOBEC.table) <- APOBEC.table$gene_name
APOBEC.table <-data.frame(t(APOBEC.table[,c(2:51)]))
sample_name<- data.frame(c("S03","S07","S08","S10","L21","L22","L27","L29","L36")) %>% `colnames<-`(c("sample"))
APOBEC.table <- APOBEC.table[grepl(paste(sample_name$sample, collapse="|"), rownames(APOBEC.table)), ]

##
APOBEC.table_raw <- row.count[rownames(row.count)%in%APOBEC$gene_id,]
APOBEC.table_raw$gene_id <- rownames(APOBEC.table_raw)
APOBEC.table_raw <- merge(APOBEC.table_raw, APOBEC, by = "gene_id")
rownames(APOBEC.table_raw) <- APOBEC.table_raw$gene_name

APOBEC.table_raw <-data.frame(t(APOBEC.table_raw[,c(2:51)]))
APOBEC.table_raw <- APOBEC.table_raw[grepl(paste(sample_name$sample, collapse="|"), rownames(APOBEC.table_raw)), ]



APOBECA.table$sample <- factor(APOBECA.table$sample, levels = APOBECA.table$sample)
ggplot(APOBECA.table, aes(x=sample, y=CHM13_G0037458))+
  geom_col(width=0.8, fill="#01B468")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(x="Sample", y="CPM")+
  ggtitle("APOBEC3A")


APOBECB.table$sample <- factor(APOBECB.table$sample, levels = APOBECB.table$sample)
ggplot(APOBECB.table, aes(x=sample, y=CHM13_G0037459))+
  geom_col(width=0.8, fill="steelblue")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(x="Sample", y="CPM")+
  ggtitle("APOBEC3B")
  
#################################################################################################
###change genelist
gene.table <- read.table("/mnt/nas2/yh/10.RNA_seq_poly_A/8.count_ucsc/gene.list2",header = F) %>%`colnames<-`(c("gene_id","gene_name"))
pd.list <- read.table("/mnt/nas2/yh/10.RNA_seq_poly_A/reference/htseq/ucsc/protein_coding.list", header = F) %>% `colnames<-`(c("gene_name"))
pd.table <- merge(gene.table,pd.list,by = "gene_name") 

cpm.xls <- data.frame(normalize_cpm.table)
cpm.xls$gene_id <- rownames(cpm.xls)
cpm.xls.table <- merge(cpm.xls,gene.table,by="gene_id")
write.csv(cpm.xls.table,"/mnt/nas2/yh/10.RNA_seq_poly_A/8.count_ucsc/cpm_20240117.xls")


cpm_pd <- data.frame(normalize_cpm.table)
  cpm_pd$gene_id <- rownames(cpm_pd)
  cpm_pd.table <- data.frame(merge(cpm_pd,pd.table,by="gene_id"))%>%dplyr::select(2:51)
  cpm_pd.table_xls <- data.frame(merge(cpm_pd,pd.table,by="gene_id")) 
    write.csv(cpm_pd.table_xls,"/mnt/nas2/yh/10.RNA_seq_poly_A/8.count_ucsc/cpm_pd_20240117.xls")
  
normalize_cpm.pd_long <-  cpm_pd.table[,grep("L", colnames(cpm_pd.table))]
  normalize_cpm.pd_long <- normalize_cpm.pd_long[, order(colnames(normalize_cpm.pd_long))]
  
normalize_cpm.pd_short <-  cpm_pd.table[,grep("S", colnames(cpm_pd.table))]
  normalize_cpm.pd_short <- normalize_cpm.pd_short[, order(colnames(normalize_cpm.pd_short))]
  
corre.table_short <- cor(normalize_cpm.pd_short ,method = 'spearman')
corre.table_long <- cor(normalize_cpm.pd_long ,method = 'spearman')
corre.table <- cor(cpm_pd.table ,method = 'spearman')
  
  round(corre.table_short, 2)
  round(corre.table_long, 2)
  round(corre.table, 2)
  
  short_melted_cormat <- melt(corre.table_short)
  long_melted_cormat <- reshape2::melt(corre.table_long)
  melted_cormat <- reshape2::melt(corre.table)
  
ggplot(data = melted_cormat, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile(color = "gray") +
    #geom_text(aes(label = round(value, digits = 2)), color = "#444444", size = 4) +
    labs(x = "", y = "", fill = "", title = "") +
    coord_fixed() +
    theme_minimal() +
    scale_fill_gradientn(
      limits = c(0.6, 1),  # Adjusted the limits to cover the full range of correlation values
      colours = c("#0066FF", "white", "#FF5511"), 
      values = scales::rescale(c(0.6, 0.85, 1))
    ) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1),
      axis.text.y = element_text(angle = 0, vjust = 1, size = 10, hjust = 1),
      legend.position =c(1.05,0.935)
    )


################################################################################################


