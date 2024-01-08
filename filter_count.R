############select gene with rowcount >5  
gene_select <- data.frame(row.count_filter)
   sample <- colnames(gene_select)
oncogene_list <- read.table("/mnt/nas2/yh/10.RNA_seq_poly_A/8.count_ucsc/oncogene.list")%>%`colnames<-`(c("gene_name"))
   gene.without_ercc <- gene.table[!gene.table$gene_id%like%"^ERCC-",]

data <- data.frame()
for (f in 1:length(sample)){
  sample_name <- sample[f]
  sample.table <- data.frame(gene_select[,colnames(gene_select)%in%sample_name])%>% 
                    `rownames<-`(rownames(gene_select)) %>%
                    `colnames<-`(sample_name)
  select_cpm <- data.frame(apply(sample.table , 2, function(x) x[x>5]))
  select_cpm$gene_id <- rownames(select_cpm)
  length(select_cpm$gene_id)
  select_gene <- merge(select_cpm , gene.without_ercc, by="gene_id")
  select_gene.list <- paste(unique(select_gene$gene_name), collapse = ',')
  gene.split <- length(unique(sort(select_gene$gene_name)))
  gene.split_name <- unique(sort(select_gene$gene_name))
  if (gene.split > 10000) {
      gene.split.list <- data.frame(c("gene"))
    for (i in 1:10) {
      num <- (gene.split%/%4000)+1
      if (num < 10 && i > num) {
          list <- c("NA")
      } else if (i == num && gene.split < (i*4000)) {
          list <- gene.split_name[(((i-1)*4000)+1):gene.split]
      } else {
          list <- gene.split_name[(((i-1)*4000)+1):(i*4000)]
      }
      coln <- paste0("col",i,collapse=',')
      table <- data.frame(paste(list, collapse = ',')) %>% `colnames<-`(c(coln))
      gene.split.list <- cbind(gene.split.list,table)
    }
      
  }

  sample_gene <- data.frame(sample_name, length(unique(sort(select_gene$gene_name))), gene.split.list)
                            
  #select_gene_onco <- merge(select_gene, oncogene_list , by="gene_name")
  #select_gene_onco.list <- paste(select_gene_onco$gene_name, collapse = ',') 
  #sample_gene_onco <- data.frame(sample_name,select_gene_onco.list,length(select_gene_onco$gene_name))
  
  data <- rbind(data,sample_gene) 
}

filt.cpm <- data
write.table(all.filt.cpm,"/mnt/nas2/yh/10.RNA_seq_poly_A/8.count_ucsc/rna_filt_gene_raw.txt",sep="\t")
