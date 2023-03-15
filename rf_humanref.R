library(randomForest)
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(dior)
setwd("/storage/chenlab/Users/qingnanl/humanre/out/humanref/")
href <- dior::read_h5("href.h5")

run_rf <- function(ref, qry, refkey, qrykey, name_label, cutoff = 0.2){
  print("preprocessing")
  qry <- FindVariableFeatures(qry, selection.method = "vst", nfeatures = 2000)
  genes.use <- intersect(VariableFeatures(qry), rownames(ref))
  print(length(genes.use))
  training.set <- c()
  test.set <- c()
  training.label <- c()
  test.label <- c()
  Idents(ref) <- refkey
  Idents(qry) <- qrykey
  print("start collect training data")
  for (i in unique(as.character(ref@meta.data[[refkey]]))){
    cells.in.clust <- WhichCells(ref,idents = i);
    n <- min(1000, round(length(cells.in.clust)*0.7))
    train.temp <- cells.in.clust[sample(length(cells.in.clust))][1:n]
    test.temp <- setdiff(cells.in.clust, train.temp)
    training.set <- c(training.set,train.temp)
    test.set<-c(test.set,test.temp)
    training.label <- c(training.label, rep(i,length(train.temp)))
    test.label <- c(test.label, rep(i, length(test.temp)));
  }
  predictor_Data <- as.matrix(ref@assays$RNA@data[genes.use,])
  tmp <- as.vector(table(training.label))
  sampsizes <- rep(min(tmp),length(tmp))
  predictor_Data <- t(scale(t(predictor_Data), center=TRUE, scale=TRUE))
  predictor_Data[is.na(predictor_Data)] <- 0
  print("start train rf")
  rf_obj <- randomForest(x=t(predictor_Data[,training.set]), y=factor(training.label), 
                          importance = TRUE, ntree = 1001, proximity=TRUE, sampsize=sampsizes, 
                          keep.inbag=TRUE, replace=FALSE)
                          
  saveRDS(rf_obj, paste0(name_label, "_rf_obj.rds"))
  print("finished train rf")
  test.predict <- predict(rf_obj,t(predictor_Data[,test.set]))
  test.predict.vote <- predict(rf_obj,t(predictor_Data[,test.set]), type = "vote")

  test.vote.cutoff<-cutoff
  test.predict.1<-c()
  for (i in 1:nrow(test.predict.vote)){
    if(max(test.predict.vote[i, ]) < test.vote.cutoff){
      test.predict.1[i]<-"unassigned"
    }else{
      test.predict.1[i]<-names(which.max(test.predict.vote[i, ]))
    }
  }
  confusion.test.1<-table(test.label, test.predict.1)
  confusion.test.1<-data.frame(confusion.test.1)
  confusion.test.1<-group_by(confusion.test.1, test.label) %>% mutate(percent = Freq/sum(Freq))
  write.csv(confusion.test.1, paste0(name_label, "_train_match.csv"))
  print("finished validation in ref dataset")
  #
  #predict
  qry.rf <- as.matrix(qry@assays$RNA@data[genes.use,])
  qry.rf <- t(scale(t(qry.rf), center=TRUE, scale=TRUE))
  qry.rf[is.na(qry.rf)] <- 0

  qry.ident <- factor(qry@meta.data[[qrykey]]) 
  names(qry.ident) <- colnames(qry)
  qry.predict.vote = predict(rf_obj,t(qry.rf), type = "vote")
  qry.vote.cutoff<-cutoff
  qry.predict.1<-c()
  for (i in 1:nrow(qry.predict.vote)){
    if(max(qry.predict.vote[i, ]) < qry.vote.cutoff){
      qry.predict.1[i]<-"unassigned"
    }else{
      qry.predict.1[i]<-names(which.max(qry.predict.vote[i, ]))
    }
  }
  confusion.qry.1 <- table(qry.ident, qry.predict.1)
  confusion.qry.1 <- data.frame(confusion.qry.1)
  confusion.qry.1 <- group_by(confusion.qry.1, qry.ident) %>% mutate(percent = Freq/sum(Freq))
  write.csv(confusion.qry.1, paste0(name_label, "_qry_match.csv"))
  print("finished")
}

setwd("/storage/chenlab/Users/qingnanl/humanre/out/humanref/")
# href <- dior::read_h5("href.h5")
Idents(href) <- "major_type"
href_BC <- subset(href, idents = "BC")
hqry_BC <- dior::read_h5("/storage/chenlab/Users/qingnanl/humanre/out/subtypes/BC_final.h5")
run_rf(ref = href_BC, qry = hqry_BC, refkey = "cell_type", qrykey = "subtype", name_label = "BC")

href_AC <- subset(href, idents = "AC")
hqry_AC <- dior::read_h5("/storage/chenlab/Users/qingnanl/humanre/out/subtypes/AC_final.h5")
run_rf(ref = href_AC, qry = hqry_AC, refkey = "cell_type", qrykey = "subtype", name_label = "AC", cutoff = 0.1)

href_HC <- subset(href, idents = "HC")
hqry_HC <- dior::read_h5("/storage/chenlab/Users/qingnanl/humanre/out/subtypes/HC_final.h5")
run_rf(ref = href_HC, qry = hqry_HC, refkey = "cell_type", qrykey = "subtype", name_label = "HC")


href_RGC <- subset(href, idents = "RGC")
Idents(href_RGC) <- "cell_type"
href_RGC <- subset(href_RGC, idents = c("MG_OFF", "MG_ON"), invert = T)
hqry_RGC <- dior::read_h5("/storage/chenlab/Users/qingnanl/humanre/out/subtypes/RGC_final.h5")
Idents(hqry_RGC) <- "subtype"
hqry_RGC <- subset(hqry_RGC, idents = c("ON_MGC", "OFF_MGC"), invert = T)

run_rf(ref = href_RGC, qry = hqry_RGC, refkey = "cell_type", qrykey = "subtype", name_label = "RGC", cutoff = 0.2)

# 

