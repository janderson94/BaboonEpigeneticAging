#!/usr/bin/env Rscript

#Read in meta file (like the one provided in this repository)
info<-read.table("./meta_info.txt",header=T)


#Plot predicted age versus chronological age for males and females
getwd()

info_n277<-read.table("~/Desktop/info_n277.txt",header=T)
info_n286<-read.table("~/Desktop/info_n286.txt",header=T)

si<-read.csv("~/Desktop/Anderson_et_al_2020_NatComm_Supplementary_Tables_1Oct20.csv")
