######################################
# Compare an imputed vcf with the original full vcf
# Record results as the number of samples with missing data
# and % of varaints removed and the % of imputed variants imputed wrongly
######################################
library(SNPRelate)

originalVCFFilename = "original.vcf"

# get file names
imputedFiles = list.files("imputed/", "*.vcf.gz", full.names = TRUE)


#create gds - only need to do this once, so check if it already exists.
gdsFilename = substr(originalVCF,1,nchar(originalVCF)-4)
gdsFilename = sprintf("%s.gds",gdsFilename)
if(!file.exists(gdsFilename))
{
  snpgdsVCF2GDS(originalVCF, gdsFilename)
}

# load geno matrix
gds = snpgdsGetGeno(gdsFilename, snpfirstdim=FALSE, with.id=TRUE)
numSNPs = length(gds$snp.id)

#initialise data frame to store results
results = data.frame(matrix(vector(), 0, 4,dimnames=list(c(), c("samples", "rate", "different.total", "different.pc"))),stringsAsFactors=F)

# check if continuing - will have to filter imputedFiles
if(file.exists("completed_vcf.csv"))
{
  alreadyDone = read.table("completed_vcf.csv")
  #filter imputedFiles
  imputedFiles = setdiff(imputedFiles, alreadyDone[,1])
  
} else
{
  #initialise missing csv
  write.table(results, "missing.csv")
}

for(i in 1:length(imputedFiles))
{
 
  # find the number of samples with missing SNPs and the number of missing SNPs
  currFileName = unlist(strsplit(imputedFiles[i], "//"))[2]
  tmp = unlist(strsplit(currFileName, "_"))
  samples = strtoi(tmp[1])
  rate = strtoi(unlist(strsplit(tmp[3], "pc"))[1])
  rate = rate/100
  missing = floor(rate * numSNPs)
  print(sprintf("samples:%i,rate:%f,missing:%d",samples,rate,missing))

  #current file to check
  #create gds
  snpgdsVCF2GDS(imputedFiles[i], "tmp/current.gds")
  print(imputedFiles[i])
  
  # load geno matrix
  gdsCheck = snpgdsGetGeno("data/tmp/current.gds", snpfirstdim=FALSE, with.id=TRUE)
  
  wrong = sum((gds$genotype[1:samples,] == gdsCheck$genotype[1:samples,])==FALSE)
  
  
  totalMissing = missing*samples
  newRow <- data.frame(samples = samples, rate = rate,  different.total = wrong, different.pc = wrong/totalMissing)
  
  results = rbind(results, newRow)
  print(results)
  
  write.table(newRow, "missing.csv", append=TRUE, col.names=FALSE, row.names=FALSE)
  write.table(imputedFiles[i], "completed_vcf.csv", append=TRUE, col.names=FALSE, row.names=FALSE)
}

write.table(results, "missing.final.csv")
write.table(imputedFiles, "completed_vcf.final.csv")

