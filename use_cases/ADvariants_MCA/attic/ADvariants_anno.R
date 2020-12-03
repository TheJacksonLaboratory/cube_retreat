library(data.table)

variants <- read.delim("ADvariants_anno_0.txt", header=TRUE, sep="\t")[,-c(2:5)]
inpath <- "/sdata/carter-lab/carter/AMPAD/synapseData/Joint/WGS"
annofiles <- list.files(path=inpath, pattern = ".annotated.txt")

anno.all <- c()
for (fl in annofiles) {
  print(fl)
  annot <- fread(file=paste(inpath, fl, sep="/"), sep="\t")
  vi <- which(annot$ID%in%variants$rsID)
  anno.all <- rbind(anno.all, annot[vi,c("ID", "CHROM","POS", "REF","ALT")])
  rm(annot)
  }

var.anno <- merge(anno.all, variants, by="ID", all=TRUE)

fwrite(var.anno, 
        file="ADvariants_anno.txt", sep="\t", 
        row.names = FALSE, quote = FALSE)
