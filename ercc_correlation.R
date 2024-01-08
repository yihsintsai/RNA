ercc_cor <- function(ercc_control.table,rna_table,condition,corr_arg){
#require file
  #ercc_control.table
  #rna_table
  #condition
  #corr_arg

lib <- c("RUVSeq","data.table","corrplot","dplyr")

for(f in lib){
  if(!require(f,character.only = TRUE)) install.packages(f)
  library(f,character.only = TRUE)
}

############## ERCC control table preparation #####################
ERCC_con_table <- data.frame(ercc_control.table) %>%
  dplyr::select(1:4) %>%
  `names<-`(c("RE-sort_ID","ERCC_ID","subgroup","Conc.Mix1"))

############## ERCC Expect molecules ##############################
# Perform calculations to ERCC copies
# transform Mix 1 concentrations= values to reflect dilutions in experiments
# 1 x 10^18 attomoles = 1 mole
# avogadros number = 6.022 x 10^23 molecules/mole
#--dilution = 1000
#--ercc_microliter  6 ul
#--rna_microgram 0.3 ug
ercc_calculations <- ERCC_con_table %>%
  mutate(Attomoles.added = (Conc.Mix1/1000)*6) %>%
  mutate(Moles.added = Attomoles.added/(10^18)) %>%
  mutate(Molecules.added = Moles.added * (6.022 * 10^23)) %>%
  #mutate(LogMolecules = log2(Molecules.added)) %>%
  dplyr::select(2,7) %>%
  rename(expect_conc=Molecules.added)

############## Filter data (> 5 reads in at least two samples for each gene) ########
len <- length(rna_table)
filter <- apply(rna_table, 1, function(x) length(x[x>5])>=2)
filtered <- rna_table[filter,]
con <- condition
x <- as.factor(con)
#grep ERCC
ERCC_counts <- filtered[row.names(filtered)%like%"^ERCC-",]
ERCC_counts$ERCC_ID <- rownames(ERCC_counts)

############## Correlation ########################################
ERCC_cor_table <- merge(ercc_calculations,ERCC_counts, by="ERCC_ID" ) %>%
  dplyr::select(3:(len+2))
arg <-c(corr_arg)
correlation <- cor(ERCC_cor_table ,method = arg)
round(correlation, 2)

col<- colorRampPalette(c("white","#FAA0A9","#FA394D"))(299)
result <- corrplot(correlation,
         method="color",
         col=col,
         number.cex = 1,
         col.lim = c(min(0.8), max(1)),
         tl.cex=1.5,
         type="full",
         order="hclust",
         addCoef.col = "black",  # Add coefficient of correlation
         tl.col="black",
         tl.srt=90,  #Text label color and rotation
         is.corr = FALSE)

}
