#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
repy = args[1] 
method=args[2]

library("ggplot2")
library("gridExtra")
rootpath=paste0("/group/qtlchenggrp/tianjing/GWAS/y",repy,"/",method)
setwd(rootpath)

realQTL=read.table(paste0("/group/qtlchenggrp/tianjing/GWAS/y",repy,"/QTL_ID.txt"))$V1
mapfile=read.csv("/group/qtlchenggrp/tianjing/GWAS/pig_map_chr1.csv")
realQTL_map = mapfile[mapfile$Column1 %in% realQTL,]

if (grepl("hmc", method, fixed=TRUE)){  #NN-Bayes
  wppa=read.table("GWAS_WPPA4Bayesnn.txt",header=T,sep=",")
} else {  #linear model
  wppa=read.table("GWAS_results_MCMC_samples_marker_effects_geno_y.txt",header=T,sep=",")
}

wppa = wppa[order(wppa$window),]

wppa$realQTL="Non_QTL"
for (i in 1:nrow(wppa)){  #each window
  window_start=wppa$start_SNP[i]
  window_end = wppa$end_SNP[i]
  cat("window ",i," start: ",window_start," ; end: ",window_end,"\n")
  for (j in 1:nrow(realQTL_map)){  #each QTL
    QTL_pos = realQTL_map$pos[j]
    if (QTL_pos>=window_start & QTL_pos<=window_end){
      cat("QTL ",realQTL_map$Column1[j]," pos ",QTL_pos," in window ",i,"\n")
      wppa$realQTL[i] = "QTL"
    }
  }
}

################################# ROC
library(pROC)

roc_y=roc(wppa$realQTL, wppa$WPPA)
auc_y=auc(roc_y)

jpeg(paste0("roc.y",repy,".",method,".jpg"))
plot(roc_y, print.auc = TRUE)
dev.off()

write.csv(auc_y,paste0("auc.y",repy,".",method,".csv"))
