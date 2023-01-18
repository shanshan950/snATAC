#! /usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
require("openxlsx")
require("reshape2")
print(args[1])
## To read in the layout file in excel format
All.barcodes<-read.xlsx(args[1],sheet=1,colNames=F,rowNames=T)
Layout.matrix<-read.xlsx(args[1],sheet=2)
Layout.matrix.mt<-melt(Layout.matrix,id.vars="X1")
BC.list<-list(
N701_N712=c("Nextera_N701","Nextera_N702","Nextera_N703","Nextera_N704","Nextera_N705","Nextera_N706","Nextera_N707","Nextera_N708","Nextera_N709","Nextera_N710","Nextera_N711","Nextera_N712"),
N713_N724=c("Nextera_N713","Nextera_N714","Nextera_N715","Nextera_N716","Nextera_N717","Nextera_N718","Nextera_N719","Nextera_N720","Nextera_N721","Nextera_N722","Nextera_N723","Nextera_N724"),
N725_N736=c("Nextera_N725","Extended_N726","Extended_N727","Extended_N728","Extended_N729","Extended_N730","Extended_N731","Extended_N732","Extended_N733","Extended_N734","Extended_N735","Extended_N736"),
N501_N508=c("Nextera_N501","Nextera_N502","Nextera_N503","Nextera_N504","Nextera_N505","Nextera_N506","Nextera_N507","Nextera_N508"),
N509_N516=c("Nextera_N509","Nextera_N510","Nextera_N511","Nextera_N512","Nextera_N513","Nextera_N514","Nextera_N515","Nextera_N516"),
N517_N524=c("Nextera_N517","Nextera_N518","Extended_N519","Extended_N520","Extended_N521","Extended_N522","Extended_N523","Extended_N524"),
N525_N532=c("Extended_N525","Extended_N526","Extended_N527","Extended_N528","Extended_N529","Extended_N530","Extended_N531","Extended_N532"),
N533_N540=c("Extended_N533","Extended_N534","Extended_N535","Extended_N536","Extended_N537","Extended_N538","Extended_N539","Extended_N540")
)
for (plate in 1:nrow(Layout.matrix.mt)){
BarcodeIndex<-All.barcodes[c(BC.list[[Layout.matrix.mt[plate,2]]],BC.list[[Layout.matrix.mt[plate,1]]]),,drop=F]
name<-paste(Layout.matrix.mt[plate,3],"BCindex",sep=".")
print(paste("Write barcode index",name))
write.table(BarcodeIndex,name,quote=F,sep="\t",col.names=F)
}
