# Comnbine and recode variants genotype
files <- list.files(pattern="ADvariants_genotypes_")
variants <- read.delim(file="ADvariants_anno.txt", colClasses = "character")
genotypesL <- c()
  
for (fn in files) {
  gfile <- read.delim(file=fn, header =FALSE, colClasses = "character")
  rsid <- gsub(pattern="ADvariants_genotypes_(rs\\d+)", "\\1", x=fn)
  g_rsid <- paste(variants$ClosetGene[variants$ID==rsID], rsID, sep="_")
  pos <- c(unique(gfile[,2]))
  genotypesL[[g_rsid]]$g_rsID <- gfile[gfile[,2]==pos,6]
  genotypesL[[g_rsid]]$patientID <- unique(as.character(gfile[,5]))
  }
 

genotypes <- as.data.frame(genotypesL[[1]])
for (l in 2:length(genotypesL)) {
  genotypes <- merge(genotypes, as.data.frame(genotypesL[[l]]), by="patientID")
  }

rownames(genotypes) <- genotypes$patientID
genotypes <- genotypes[,-1]

genotypes <- ifelse(genotypes=="0/0", "0", ifelse(genotypes=="0/1", "1", "2"))

## Combine APOE variants to one genotype APOE_e4: 0 for e2/e2; 1 for hets; 2 for e4/e4

###

write.table(genotypes, file="ADvariants_genotypes_all.csv", quote=F, sep=",")
