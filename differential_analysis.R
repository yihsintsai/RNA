
title: "differential analysis for long survival"
author: "yh"
date: "2024-02-27"
output: html_document

library("edgeR")
library("RUVSeq")
library("data.table")
library("RColorBrewer")
library("dplyr")
library("ggcorrplot")
library("DESeq2")
library(devtools)
library(ggvenn)
library("pdist")
library(enrichR)
library("stringr")
source("/mnt/nas2/yh/10.RNA_seq_poly_A/3.R/ercc_cor.R")
setwd("/mnt/nas2/yh/10.RNA_seq_poly_A/8.count_ucsc/")

## data prepare
all_countData <- counts(ercc_data)
  long_count <- all_countData[,grep("L21|L23|L24|L26|L30|L31|L32|L35",colnames(all_countData))]
all_colData <- pData(ercc_data)
  long_colData <- data.frame(all_colData[grep("L21|L23|L24|L26|L30|L31|L32|L35",rownames(all_colData)),])
    long_colData$condition <- c("OP1","OP2","OP1","OP2","OP1","OP2","OP3","OP1","OP2","OP3",
                            "OP1","OP2","OP3","OP1","OP2","OP1","OP2","OP1","OP2")


long_design <- model.matrix(~ W_1 + W_2 + W_3 + condition, data = long_colData)

  Long_survival_dds <- DESeqDataSetFromMatrix(countData = long_count ,
                                     colData = long_colData,
                                     design = ~ W_1 + W_2 + W_3 + condition )
## filter the sum of rowcount < 10
keep_gene <- rowSums(counts(Long_survival_dds)) >= 10
   Long_survival_dds <- Long_survival_dds[keep_gene,]
    Long_survival_dds$condition <- factor(Long_survival_dds$condition, levels = c("OP1","OP2","OP3"))

## calculate normalize count
Long_survival_dds <- estimateSizeFactors(Long_survival_dds)
  sizeFactors(Long_survival_dds)
  Long_survival_nor_counts <- data.frame(counts(Long_survival_dds, normalized=TRUE))
  Long_survival_nor_counts$gene_id <- rownames(Long_survival_nor_counts)


##cofficient of normalize 
##deseq_cofe <-list(environment(L1v2_dds@dispersionFunction)[["fit"]][["coefficients"]])


## DEG analysis
Long_survival_ds <- DESeq(Long_survival_dds)
  head(Long_survival_ds)

## construct condition list
resultsNames(Long_survival_ds)
Long_survival_res <- results(Long_survival_ds)
  condition.list <- data.frame("OP1vs.OP2"= c("OP1","OP2"),
                               "OP1vs.OP3"= c("OP1","OP3"),
                               "OP2vs.OP3"= c("OP2","OP3"))

## result table 
res_data <- data.frame(Long_survival_res@rownames) %>%
            `colnames<-`(c("gene_id"))

for (f in 1:3){
  ##condition extract
  comparison_1 <- condition.list[f][1,]
  comparison_2 <- condition.list[f][2,]
    condition <- paste(comparison_1,"vs.",comparison_2,sep="")
  
  ##differential analysis result
  res <- results(Long_survival_ds,
                 alpha=0.05,contrast=c("condition",comparison_2,comparison_1 ),
                 independentFiltering = FALSE,
                 cooksCutoff=FALSE)
  
  ##table build
  res_table <- data.frame(res)
  column.name <- colnames(res_table)
    column.name <- paste(condition,column.nmae,sep="_")
    colnames(res_table) <- c(column.name)
  res_table$gene_id <- rownames(res)
  rownames(res_table) <-  NULL 
  res_data <- merge(res_data,res_table,by = "gene_id")
}
rownames(res_data) <- res_data$gene_id

merge_table <- merge(Long_survival_nor_counts,res_data,by="gene_id") 

head(merge_table) 


## pvalue_histrogram
for (f in 1:3){
  ##condition extract
  comparison <- colnames(condition.list)[f]
    condition_title <- paste( comparison,"pvalue",sep="_")
    condition_pvalue <- res_data[,grepl(condition_title,colnames(res_data))]
    pvalue_histrogram <- hist(condition_pvalue,
                               breaks=100, 
                              col="grey" , 
                              main= paste("Histrogram of Long survival",comparison,sep=" "))
  
  
}


## merge gene location

gene_location <- read.csv("/mnt/nas2/elis/ref/chm13.draft_v2.0.gene_annotation.sorted.filtered.vfilter.csv")
  gene_location.table <- merge(gene_location,merge_table,by = "gene_id" )
  head(gene_location.table)

write.csv(gene_location.table,"/mnt/nas2/yh/10.RNA_seq_poly_A/8.count_ucsc/Long_survival_DE_analysis_20240229.xlsx")

  ##check loose gene
  loss_gene <- merge_table[!merge_table$gene_id%in%gene_location.table$gene_id,]
    print(nrow(loss_gene))
    head(loss_gene$gene_id)


## vocano plot

cut_off <- 0.01
significant.table <- data.frame(row.names =c("Up-regulated (fdr <0.01 & lof2FC >0)",
                                             "Down-regulated (fdr <0.01 & lof2FC >0)",
                                             "No-regulated") )
gene.list <- list()
for (f in 1:3){
  comparison <- colnames(condition.list)[f]
  
  data <- gene_location.table[,grepl(comparison,colnames(gene_location.table))]%>%`rownames<-`(gene_location.table$gene_id)
    colnames(data) <- gsub(paste(comparison,"_",sep=""),"",colnames(data))
    data$expressed <- "NO"
    data$expressed[data$log2FoldChange > 0 & data$padj < cut_off] <- "UP"
    data$expressed[data$log2FoldChange < 0 & data$padj < cut_off] <- "DOWN"

    significant_gene_up <- rownames(data[(data$expressed%in%"UP"),])
    significant_gene_down <- rownames(data[(data$expressed%in%"DOWN"),])
    significant_gene <- list(gene_list_UP = c(significant_gene_up),gene_list_DOWN=c(significant_gene_down))
    gene.list <- c(gene.list,significant_gene)
    names(gene.list) <- gsub("gene_list", comparison, names(gene.list))
        
    count <- data.frame(c(NROW(grep("UP",data$expressed)),
                          NROW(grep("DOWN",data$expressed)),
                          NROW(grep("NO",data$expressed)))) %>%
                          `rownames<-`(c("Up-regulated (fdr <0.01 & lof2FC >0)",
                                         "Down-regulated (fdr <0.01 & lof2FC >0)",
                                         "No-regulated")) %>%
                          `colnames<-`(comparison)
  
    significant.table <- cbind(significant.table,count) 
    

  vocano_fig <- ggplot(data = data, aes(x = log2FoldChange, y = (-log10(padj)) , col = expressed)) +
                        geom_vline(xintercept = c(0), col = "gray", linetype = 'dashed') +
                        geom_hline(yintercept = -log10(0.01), col = "gray", linetype = 'dashed') + 
                        geom_point() + 
                        theme_minimal() +
                        labs(x="Log2FoldChange", y="-Log10(padj)", title=paste("Long survival", comparison , sep =" ")) +
                        theme_set(theme_classic(base_size = 20) +
                        theme(
                              axis.title.y = element_text (face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
                              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
                              plot.title = element_text(hjust = 0.5,face = "bold"))
                                 ) +
                        scale_color_manual(values = c("blue","gray", "#bb0c00"), 
                                           labels = c("Down-regulated","No-regulation", "Up-regulated" ))
  plot(vocano_fig)

}
print(significant.table)
summary(gene.list)


## gene of comparison intersect

inters <-gene.list
myCol <- brewer.pal(12, "Paired")[c(1,5,3,7)]


ggvenn(
  inters, 
  fill_color = myCol ,
  stroke_size = 0.3, set_name_size = 4
)

##split up&down
inters <-c(
  OP1vs.OP2 = reshape2::melt(gene.list[grepl("OP1vs.OP2",names(gene.list))])[1],
  OP1vs.OP3 = reshape2::melt(gene.list[grepl("OP1vs.OP3",names(gene.list))])[1],
  OP2vs.OP3 = reshape2::melt(gene.list[grepl("OP2vs.OP3",names(gene.list))])[1]
)
ggvenn(
  inters, 
  fill_color = myCol ,
  stroke_size = 0.3, set_name_size = 4
)


## cat union gene
    
col <- list(sample_id = c(OP1= "#02DF82", 
                          OP2="#66B3FF",
                          OP3="#FFA042"))
clust_col <- colorRampPalette(c("blue","white","#FA394D"))(299)

tree.table <- list()

union.gene <- unique(reshape2::melt(inters)[1])
  write.csv(union.gene, "/mnt/nas2/yh/10.RNA_seq_poly_A/8.count_ucsc/union.gene.csv")
union.gene.table <- Long_survival_nor_counts[Long_survival_nor_counts$gene_id%in%union.gene$value,]
  rownames(union.gene.table) <- union.gene.table$gene_id

clust_table<-c("union.gene.sort.table","union.gene.unsort.table ")
for (f in 1:length(clust_table)){
  type <- clust_table[f]
  if ( type == "union.gene.sort.table") {
    clust.table <- union.gene.table[,order(gsub("L.*_","",colnames(union.gene.table)))][c(-1)]
    clust.tree.table <- clust.table
  } else {
    clust.table <- union.gene.table[,order(colnames(union.gene.table))][c(-1)]
  }
   
    annotation_col = data.frame(sample_id=factor(gsub("L.*_","",colnames(clust.table))))%>%
                                `rownames<-`(colnames(clust.table))
    
  if ( type == "union.gene.sort.table") {
      
     cluster_tree <- pheatmap(clust.table,
                              col= clust_col,
                              clustering_method = "ward.D2",
                              cluster_row =T,
                              cluster_col =F,
                              scale = "row", 
                              border = T ,
                              border_color = "black", 
                              gaps_col =  c(2,4,7,10,13,15,17),
                              show_rownames = F,
                              width = 13,
                              hight = 13
                              )
     tree.table <- c(tree.table ,list(union.gene.sort=cluster_tree[["tree_row"]]))
      
    } else {
    cluster_tree <- pheatmap(clust.table,
                             col= clust_col,
                             clustering_method = "ward.D",
                             cluster_row =T,
                             annotation_col = annotation_col,
                             annotation_colors = col,
                             cluster_col =F,
                             scale = "row",
                             border = T ,
                             border_color = "black",
                             show_rownames = F,
                             width = 13,
                             hight = 13)
    tree.table <- c(tree.table ,list(union.gene.unsort=cluster_tree[["tree_row"]]))
  }
    
}


#names <- rownames(geneExp_matrix)


## group from clust tree

tree <- tree.table[["union.gene.sort"]]
  clust_tree <- data.frame(cutree(tree, k=4) ) %>% `colnames<-`(c("group"))
    clust_tree$gene_name <- rownames(clust_tree)
    
clust.data <- list()
OP.group <- c("OP1","OP2","OP3")
plot.list <- list()
for (f in 1:4){
  name <- paste("group",f, sep = "_")
  clust.name <- clust_tree[(clust_tree$group %in% f),] %>% dplyr::select(2)
    print(paste(name,(length(clust.name$gene_name)),sep=":"))
    
  data <- data.frame(clust.tree.table[rownames(clust.tree.table)%in%clust.name$gene_name,])
  
  group.data <- list(data = data, count=length(clust.name$gene_name))
  
  table <- data.frame(rownames(data))%>%`colnames<-`(c("gene_id"))
  for (i in 1:length(OP.group)) {
    OP = OP.group[i]
    group <- data[,grep( OP , (colnames(data))) ]    
      group$mean <- rowMeans(group)
      group$gene_id<- rownames(group)
      group <-  group[,!grepl( OP , (colnames(group)))]
    table <- merge(table,group,by="gene_id")
      names(table) <- gsub("mean",OP, names(table)) 
  }
  group.data <- c(group.data,list(mean=table))
  clust.data <- c(clust.data,list(data=group.data))
  names(clust.data) <- gsub("data",name, names(clust.data))

##plot table  
  rownames(table)<-table$gene_id
  table <- data.frame(scale(t(table[,grepl("OP",colnames(table))])))
    table$sample <- rownames(table)
    table <- table %>% gather(key="CHM", value='value', -sample)
  plot.list <- c(plot.list,list(data = table))
    names(plot.list) <- gsub("data",name, names(plot.list))
}


## plot mean figure

for (f in 1:length(plot.list)){
  group <- plot.list[[f]]
  
  color=rainbow(n = 4,s=.8, v=.9)
    color=color[f]
  theme_set(theme_minimal())
  
  mean_fig <- ggplot(group, aes(x = sample, y = value, group=CHM))+
                    geom_line(color=color)+
                    geom_point()+
                    #geom_smooth(method = "lm", se = FALSE)+
                    scale_x_discrete()+
                    scale_color_brewer()+
                    theme_minimal()+
                    theme_bw() +  
                    theme(legend.position = "none", 
                          axis.text.x = element_text(angle = 90 , vjust = 0.4)) +
                    labs(title = paste(names(plot.list)[f],length(unique(group$CHM)),sep=" : "))
                    
 
 
 plot(mean_fig) 
}

## GO database


db_library <- c("GO_Biological_Process_2023","GO_Cellular_Component_2023","GO_Molecular_Function_2023")

websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes   
  }
if (websiteLive) dbs <- listEnrichrDbs()
if (websiteLive) head(dbs)

db_library <- c("GO_Biological_Process_2023","GO_Cellular_Component_2023","GO_Molecular_Function_2023")


  
## GO plot
library("stringr")    
gene_ref <- read.table("/mnt/nas2/yh/10.RNA_seq_poly_A/reference/go_reference/entrez_ID_722_genes.txt", header =T,sep="\t" )

enrich.table <- list()

for (i in 1:length(plot.list)){
  gene <- clust.data[[i]][['mean']][["gene_id"]]
  gene_name <- gene_ref[gene_ref$gene_ID %in% gene,] %>% select(contains("gene_name"))
  go.table <- list()                                                        
  for (f in 1:length(db_library)){
    dbs <- db_library[f]
    enriched <- enrichr(gene_name$gene_name, dbs)
    table <- enriched[[1]]
    go.table <- c(go.table,list(data=table))
    names(go.table) <- gsub("data",dbs,names(go.table))
      
  }
  enrich.table <- c(enrich.table,list(data=go.table))
     names(enrich.table) <- gsub("data",names(clust.data)[i],names(enrich.table))
  
}

for (i in 1:length(enrich.table)){
    group <- enrich.table[[i]]
    go_plot <-data.frame ()
  for (f in 1:length(group)) {
    table <- data.frame(group[f])
      colnames(table) <- c(gsub(".*_2023.","" , colnames(table)))
      table$mapped_gene <- gsub("/.*", "", table$Overlap)
      table$all_gene <- gsub(".*/", "", table$Overlap)
      table$ratio <- as.numeric(table$mapped_gene)/as.numeric(table$all_gene)
      table <- table[order(table$P.value,decreasing = F),]
    if ( as.numeric(nrow(table)) < 15 ){
            plot_table <- table
          } else {
            plot_table <- table[c(1:15),]
          }
    plot_table$type <- c(names(group)[f])
    go_plot <-dplyr::bind_rows(go_plot, plot_table)
    
  }

    
    
    jpeg(file = (paste( "/mnt/nas2/yh/",names(clust.data)[i],"_GO.jpg",sep="" )), width = 1500, height = 1000)

       
      
      plot.fig <- ggplot(go_plot, aes(x = -log(P.value) , y = factor(Term , level=c(rev(Term))) , fill = type)) + 
                    geom_bar(stat = "identity", width = 0.7 ) +
                    scale_fill_manual(values=c("#02DF82", "#66B3FF","#FFA042"))+
                    scale_y_discrete(name = "",labels = function(x) str_wrap(x, width = 100)) +
                    scale_x_continuous(name = "-log(p-value)")+
                                      #limits = c(0,max(-log(go_plot$P.value))+5) , 
                                      #breaks = c(0:(max(-log(go_plot$P.value))+5)), minor_breaks = NULL) +
                    labs(title = (names(clust.data)[i])) +
                    geom_text(aes(label=Overlap), vjust=0.5, hjust=-0.1, size = 4 , fontface = "plain") +
                    theme_bw() + 
                    theme(axis.title = element_text(size=10,face = "plain"),
                          axis.text = element_text(color = "black", size = 10, angle = 0, hjust = 0.1, vjust = 0.5, face = "plain"),
                          panel.background = element_rect(),
                          panel.grid.major = element_line())
                     
      plot(plot.fig)
      dev.off()
}
    

