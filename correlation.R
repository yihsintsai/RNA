###change genelist
gene.table <- read.table("/mnt/nas2/yh/10.RNA_seq_poly_A/8.count_ucsc/gene.list2",header = F) %>%`colnames<-`(c("gene_id","gene_name"))
pd.list <- read.table("/mnt/nas2/yh/10.RNA_seq_poly_A/reference/htseq/ucsc/protein_coding.list", header = F) %>% `colnames<-`(c("gene_name"))
pd.table <- merge(gene.table,pd.list,by = "gene_name") 

cpm_pd <- data.frame(normalize_cpm.table)
  cpm_pd$gene_id <- rownames(normalize_cpm.table)
  cpm_pd.table <- data.frame(merge(cpm_pd,pd.table,by="gene_id"))
    %>%dplyr::select(2:51)
  cpm_pd.table_xls <- data.frame(merge(cpm_pd, pd.table, by="gene_id")) 
    write.csv(cpm_pd.table_xls,"/mnt/nas2/yh/10.RNA_seq_poly_A/8.count_ucsc/cpm_pd.xls")
  
normalize_cpm.pd_long <-  cpm_pd.table[,grep("L", colnames(cpm_pd.table))]
  normalize_cpm.pd_long <- normalize_cpm.pd_long[, order(colnames(normalize_cpm.pd_long))]
  
normalize_cpm.pd_short <-  cpm_pd.table[,grep("S", colnames(cpm_pd.table))]
  normalize_cpm.pd_short <- normalize_cpm.pd_short[, order(colnames(normalize_cpm.pd_short))]

corre.table_short <- cor(normalize_cpm.pd_short ,method = 'spearman')
   round(corre.table_short, 2)
corre.table_long <- cor(normalize_cpm.pd_long ,method = 'spearman')
  round(corre.table_long, 2)
corre.table <- cor(cpm_pd.table ,method = 'spearman')
  round(corre.table, 2)

short <- ggcorrplot(corre.table_short,lab = TRUE) +
  scale_fill_gradient2(low = "white", high = "red",midpoint = 0.65, breaks=c(0.7, 1), limit=c(0.7, 1))+
  theme(
    text = element_text(size = 12 , color = "black"), 
    axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1), 
    legend.position = c(1.03,0.95),
    legend.title = element_blank()
  )
plot(short)

long <- ggcorrplot(corre.table_long,lab = TRUE) +
  scale_fill_gradient2(low = "white", high = "red",midpoint = 0.65, breaks=c(0.6, 1), limit=c(0.6, 1))+
  theme(
    text = element_text(size = 12 , color = "black"), 
    axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),  
    legend.position = c(1.03,0.95),
    legend.title = element_blank()
  )
plot(long)

all <- ggcorrplot(corre.table,lab = TRUE) +
  scale_fill_gradient2(low = "white", high = "red",midpoint = 0.65, breaks=c(0.6, 1), limit=c(0.6, 1))+
  theme(
    text = element_text(size = 20 , color = "black"),  
    axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
    legend.position = c(1.03,0.97),
    legend.title = element_blank()
  )
plot(all)

