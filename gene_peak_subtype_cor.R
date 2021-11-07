
#gene_table is a data frame with each row as a gene and each column as a cell type (averaged from single cell by cell type annotation)
gene_table<-read.csv("AverageExpressionBySubtype.csv", row.names = 1)
df_gene<-data.frame(gene = rownames(gene_table), proportion = rowSums(gene_table != 0)/ncol(gene_table))#remove super low genes
df_gene$gene<-as.character(df_gene$gene)
rownames(df_gene)<-df_gene$gene
df_gene<-df_gene[df_gene$proportion >0.8, ]
gene_table<-gene_table[df_gene$gene, ]

#same thing, each row is a peak and each column is a cell type
peak_table<-read.csv("AveragePeakBySubtype.csv", row.names = 1)
common_col<-intersect(colnames(gene_table), colnames(peak_table))
gene_table<-gene_table[, common_col]
peak_table<-peak_table[, common_col]
df_peak<-data.frame(peak = rownames(peak_table), proportion = rowSums(peak_table != 0)/ncol(peak_table))
df_peak$peak<-as.character(df_peak$peak)
rownames(df_peak)<-df_peak$peak
df_peak<-df_peak[df_peak$proportion >0.8, ]
peak_table<-peak_table[df_peak$peak, ]

peak_table$chr<-gsub("\\:.*","",rownames(peak_table))
peak_table$start<-gsub(".*:", "", gsub("\\-.*","", rownames(peak_table)))
peak_table$end<-gsub(".*-", "",rownames(peak_table))
peak_table$start<-as.numeric(peak_table$start)
peak_table$end<-as.numeric(peak_table$end)

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(org.Hs.eg.db)
library(clusterProfiler)
gene_info<-genes(txdb)
gene_info<-as.data.frame(gene_info)
gene_info_conv<-bitr(rownames(gene_info), fromType="ENTREZID", toType=c("SYMBOL"), OrgDb="org.Hs.eg.db")
#substitute with gene names
gene_info<-gene_info[gene_info_conv$ENTREZID, ]
rownames(gene_info)<-gene_info_conv$SYMBOL
gene_info$gene_name<-rownames(gene_info)
gene_info$tss<-ifelse(gene_info$strand == "+", gene_info$start, gene_info$end)
saveRDS(gene_info, "gene_info_tss.rds")#save
common_genes<-intersect(rownames(gene_table), rownames(gene_info))#13888
write.csv(gene_table, "gene_table.csv")
write.csv(peak_table, "peak_table.csv")
#run correlation: spearman
list<-list()
k<-1
for (i in 1:length(common_genes)){
  tryCatch({
  gene_name<-as.character(common_genes[i])
  chr<-as.character(gene_info[gene_name, 1])
  tss<-gene_info[gene_name, 8]
  peak_subset<-peak_table[peak_table$chr == chr, ]
  peak_subset<-peak_subset[abs((peak_subset$start + peak_subset$end)/2 - tss) < 250000, ]
  peaks.id = rownames(peak_subset)
  models = lapply(peaks.id, function(t){
    result<-cor.test(as.numeric(peak_subset[t, 1:51]), method = "spearman", as.numeric(gene_table[gene_name, ])) #here, use 1:51 because we only compute the correlation for the first 51 column (52-54 will be chr, start, end)
    data.frame(peak = t, gene = gene_name, result[c("estimate","p.value")], tss = tss,
               stringsAsFactors=FALSE)
  		})
  models <- do.call(rbind, models);
  print(head(models))
  list[[k]]<-models
  print(i)
  print(k)
  k <- k+1
  }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
correlation_table<-do.call(rbind, list)
write.csv(correlation_table, "cor_table_sub_type_spearman.csv")
#do FDR
data<-read.csv("cor_table_sub_type_spearman.csv", row.names = 1)
head(data)
data<-data[complete.cases(data), ]
rownames(data)<-1:nrow(data)
data$chr<-gsub("\\:.*","",data$peak)
data$start<-gsub(".*:", "", gsub("\\-.*","", data$peak))
data$end<-gsub(".*-", "",data$peak)
data$start<-as.numeric(data$start)
data$end<-as.numeric(data$end)
data$distance<-abs((data$start + data$end)/2 - data$tss)
fdr<-lapply(1:nrow(data), function(t){
  gene<-as.character(data[t, "gene"])
  print(gene)
  count<-length(which(data$gene == gene))
  p.adjust(data[t, "p.value"], method = 'fdr', n = count)
})
fdr_df<-do.call(rbind, fdr)
data$fdr<-fdr_df[, 1]
write.csv(data, "cor_table_sub_type_spearman_fdr.csv")

#permutation#if you don't want to run the permutation, you can skip it
#shuffle ct
list<-list()
k<-1
for (i in 1:length(common_genes)){
  tryCatch({
  gene_name<-as.character(common_genes[i])
  chr<-as.character(gene_info[gene_name, 1])
  tss<-gene_info[gene_name, 8]
  peak_subset<-peak_table[peak_table$chr == chr, ]
  peak_subset<-peak_subset[abs((peak_subset$start + peak_subset$end)/2 - tss) < 250000, ]
  peaks.id = rownames(peak_subset)
  models = lapply(peaks.id, function(t){
    result<-cor.test(as.numeric(peak_subset[t, sample(1:51)]), method = "spearman", as.numeric(gene_table[gene_name, sample(1:51)]))
    data.frame(peak = t, gene = gene_name, result[c("estimate","p.value")], tss = tss,
               stringsAsFactors=FALSE)
  		})
  models <- do.call(rbind, models);
  print(head(models))
  list[[k]]<-models
  print(i)
  print(k)
  k <- k+1
  }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
correlation_table_permute_gene<-do.call(rbind, list)
write.csv(correlation_table_permute_gene, "cor_table_sub_type_spearman_shuffle_ct.csv")
#
data<-read.csv("cor_table_sub_type_spearman_shuffle_ct.csv", row.names = 1)
head(data)
data<-data[complete.cases(data), ]
rownames(data)<-1:nrow(data)
data$chr<-gsub("\\:.*","",data$peak)
data$start<-gsub(".*:", "", gsub("\\-.*","", data$peak))
data$end<-gsub(".*-", "",data$peak)
data$start<-as.numeric(data$start)
data$end<-as.numeric(data$end)
data$distance<-abs((data$start + data$end)/2 - data$tss)
fdr<-lapply(1:nrow(data), function(t){
  gene<-as.character(data[t, "gene"])
  print(gene)
  count<-length(which(data$gene == gene))
  p.adjust(data[t, "p.value"], method = 'fdr', n = count)
})
fdr_df<-do.call(rbind, fdr)
data$fdr<-fdr_df[, 1]
write.csv(data, "cor_table_sub_type_spearman_shuffle_ct_fdr.csv")

#
#shuffle regions
#shuffle cell types
#permutation
list<-list()
k<-1
for (i in 1:length(common_genes)){
  tryCatch({
  gene_name<-as.character(common_genes[i])
  chr<-as.character(gene_info[gene_name, 1])
  tss<-gene_info[gene_name, 8]
  peak_subset<-peak_table[peak_table$chr == chr, ]
  peak_subset<-peak_subset[abs((peak_subset$start + peak_subset$end)/2 - tss) < 250000, ]
  peak_subset<-peak_table[!(rownames(peak_table) %in% rownames(peak_subset)), ]
  x <- sample(1:nrow(peak_subset), 50)
  peak_subset<-peak_subset[x, ]
  peaks.id = rownames(peak_subset)
  models = lapply(peaks.id, function(t){
    result<-cor.test(as.numeric(peak_subset[t, 1:51]), method = "spearman", as.numeric(gene_table[gene_name, 1:51]))
    data.frame(peak = t, gene = gene_name, result[c("estimate","p.value")], tss = tss,
               stringsAsFactors=FALSE)
  		})
  models <- do.call(rbind, models);
  print(head(models))
  list[[k]]<-models
  print(i)
  print(k)
  k <- k+1
  }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
correlation_table_permute_gene<-do.call(rbind, list)
write.csv(correlation_table_permute_gene, "cor_table_sub_type_spearman_shuffle_region_2.csv")
#
data<-read.csv("cor_table_sub_type_spearman_shuffle_region_2.csv", row.names = 1)
head(data)
data<-data[complete.cases(data), ]
rownames(data)<-1:nrow(data)
data$chr<-gsub("\\:.*","",data$peak)
data$start<-gsub(".*:", "", gsub("\\-.*","", data$peak))
data$end<-gsub(".*-", "",data$peak)
data$start<-as.numeric(data$start)
data$end<-as.numeric(data$end)
data$distance<-abs((data$start + data$end)/2 - data$tss)
fdr<-lapply(1:nrow(data), function(t){
  gene<-as.character(data[t, "gene"])
  print(gene)
  count<-length(which(data$gene == gene))
  p.adjust(data[t, "p.value"], method = 'fdr', n = count)
})
fdr_df<-do.call(rbind, fdr)
data$fdr<-fdr_df[, 1]
write.csv(data, "cor_table_sub_type_spearman_shuffle_region_fdr_2.csv")


##```````````````````````````````````````````````````````````````````````````````````````````````````````
#run correlation: pearson
list<-list()
k<-1
for (i in 1:length(common_genes)){
  tryCatch({
  gene_name<-as.character(common_genes[i])
  chr<-as.character(gene_info[gene_name, 1])
  tss<-gene_info[gene_name, 8]
  peak_subset<-peak_table[peak_table$chr == chr, ]
  peak_subset<-peak_subset[abs((peak_subset$start + peak_subset$end)/2 - tss) < 250000, ]
  peaks.id = rownames(peak_subset)
  models = lapply(peaks.id, function(t){
    result<-cor.test(as.numeric(peak_subset[t, 1:51]), as.numeric(gene_table[gene_name, ]))
    data.frame(peak = t, gene = gene_name, result[c("estimate","p.value")], tss = tss,
               stringsAsFactors=FALSE)
  		})
  models <- do.call(rbind, models);
  print(head(models))
  list[[k]]<-models
  print(i)
  print(k)
  k <- k+1
  }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
correlation_table<-do.call(rbind, list)
write.csv(correlation_table, "cor_table_sub_type_pearson.csv")
#do FDR
data<-read.csv("cor_table_sub_type_pearson.csv", row.names = 1)
head(data)
data<-data[complete.cases(data), ]
rownames(data)<-1:nrow(data)
data$chr<-gsub("\\:.*","",data$peak)
data$start<-gsub(".*:", "", gsub("\\-.*","", data$peak))
data$end<-gsub(".*-", "",data$peak)
data$start<-as.numeric(data$start)
data$end<-as.numeric(data$end)
data$distance<-abs((data$start + data$end)/2 - data$tss)
fdr<-lapply(1:nrow(data), function(t){
  gene<-as.character(data[t, "gene"])
  print(gene)
  count<-length(which(data$gene == gene))
  p.adjust(data[t, "p.value"], method = 'fdr', n = count)
})
fdr_df<-do.call(rbind, fdr)
data$fdr<-fdr_df[, 1]
write.csv(data, "cor_table_sub_type_pearson_fdr.csv")

#permutation
#shuffle ct
list<-list()
k<-1
for (i in 1:length(common_genes)){
  tryCatch({
  gene_name<-as.character(common_genes[i])
  chr<-as.character(gene_info[gene_name, 1])
  tss<-gene_info[gene_name, 8]
  peak_subset<-peak_table[peak_table$chr == chr, ]
  peak_subset<-peak_subset[abs((peak_subset$start + peak_subset$end)/2 - tss) < 250000, ]
  peaks.id = rownames(peak_subset)
  models = lapply(peaks.id, function(t){
    result<-cor.test(as.numeric(peak_subset[t, sample(1:51)]), as.numeric(gene_table[gene_name, sample(1:51)]))
    data.frame(peak = t, gene = gene_name, result[c("estimate","p.value")], tss = tss,
               stringsAsFactors=FALSE)
  		})
  models <- do.call(rbind, models);
  print(head(models))
  list[[k]]<-models
  print(i)
  print(k)
  k <- k+1
  }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
correlation_table_permute_gene<-do.call(rbind, list)
write.csv(correlation_table_permute_gene, "cor_table_sub_type_pearson_shuffle_ct.csv")
#
data<-read.csv("cor_table_sub_type_pearson_shuffle_ct.csv", row.names = 1)
head(data)
data<-data[complete.cases(data), ]
rownames(data)<-1:nrow(data)
data$chr<-gsub("\\:.*","",data$peak)
data$start<-gsub(".*:", "", gsub("\\-.*","", data$peak))
data$end<-gsub(".*-", "",data$peak)
data$start<-as.numeric(data$start)
data$end<-as.numeric(data$end)
data$distance<-abs((data$start + data$end)/2 - data$tss)
fdr<-lapply(1:nrow(data), function(t){
  gene<-as.character(data[t, "gene"])
  print(gene)
  count<-length(which(data$gene == gene))
  p.adjust(data[t, "p.value"], method = 'fdr', n = count)
})
fdr_df<-do.call(rbind, fdr)
data$fdr<-fdr_df[, 1]
write.csv(data, "cor_table_sub_type_pearson_shuffle_ct_fdr.csv")

#
#shuffle regions
#shuffle cell types
#permutation
list<-list()
k<-1
for (i in 1:length(common_genes)){
  tryCatch({
  gene_name<-as.character(common_genes[i])
  chr<-as.character(gene_info[gene_name, 1])
  tss<-gene_info[gene_name, 8]
  peak_subset<-peak_table[peak_table$chr == chr, ]
  peak_subset<-peak_subset[abs((peak_subset$start + peak_subset$end)/2 - tss) < 250000, ]
  peak_subset<-peak_table[!(rownames(peak_table) %in% rownames(peak_subset)), ]
  x <- sample(1:nrow(peak_subset), 50)
  peak_subset<-peak_subset[x, ]
  peaks.id = rownames(peak_subset)
  models = lapply(peaks.id, function(t){
    result<-cor.test(as.numeric(peak_subset[t, 1:51]), as.numeric(gene_table[gene_name, 1:51]))
    data.frame(peak = t, gene = gene_name, result[c("estimate","p.value")], tss = tss,
               stringsAsFactors=FALSE)
  		})
  models <- do.call(rbind, models);
  print(head(models))
  list[[k]]<-models
  print(i)
  print(k)
  k <- k+1
  }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
correlation_table_permute_gene<-do.call(rbind, list)
write.csv(correlation_table_permute_gene, "cor_table_sub_type_pearson_shuffle_region_2.csv")
#
data<-read.csv("cor_table_sub_type_pearson_shuffle_region_2.csv", row.names = 1)
head(data)
data<-data[complete.cases(data), ]
rownames(data)<-1:nrow(data)
data$chr<-gsub("\\:.*","",data$peak)
data$start<-gsub(".*:", "", gsub("\\-.*","", data$peak))
data$end<-gsub(".*-", "",data$peak)
data$start<-as.numeric(data$start)
data$end<-as.numeric(data$end)
data$distance<-abs((data$start + data$end)/2 - data$tss)
fdr<-lapply(1:nrow(data), function(t){
  gene<-as.character(data[t, "gene"])
  print(gene)
  count<-length(which(data$gene == gene))
  p.adjust(data[t, "p.value"], method = 'fdr', n = count)
})
fdr_df<-do.call(rbind, fdr)
data$fdr<-fdr_df[, 1]
write.csv(data, "cor_table_sub_type_pearson_shuffle_region_fdr_2.csv")
