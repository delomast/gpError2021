# error rate simulation for SNP and microhaps for grandparentage with smaller panels
library(tidyverse)
library(gRandma)
library(EFGLmh)
library(parallel)

numCores <- detectCores()
if(Sys.info()["sysname"] == "Windows") numCores <- 1

### This first section is a round-about way of creating a template gRandma object
###   it is a relic of other explorations that was convenient to just copy and paste here

# some data
oxbow <- readInData("Oxbow_mh.txt")

# select most informative SNP from each locus
snpPos <- c()
for(l in getLoci(oxbow)){
	g <- oxbow$genotypes %>% select(contains(l))
	g <- unlist(g)
	g <- g[!is.na(g)]
	maf <- rep(0, nchar(g[1]))
	for(i in 1:nchar(g[1])){
		snps <- substr(g,i,i)
		if(length(unique(snps)) != 2) next # only biallelic
		maf[i] <- min(table(snps) / length(snps))
	}
	snpPos <- c(snpPos, which.max(maf))
}
names(snpPos) <- getLoci(oxbow)
oxbowSNP <- oxbow
tg <- oxbowSNP$genotypes
for(l in getLoci(oxbow)){
	tg[[paste0(l, ".A1")]] <- substr(tg[[paste0(l, ".A1")]], snpPos[l], snpPos[l])
	tg[[paste0(l, ".A2")]] <- substr(tg[[paste0(l, ".A2")]], snpPos[l], snpPos[l])
}
oxbowSNP$genotypes <- tg
oxbowSNP <- oxbowSNP %>% removeLoci("Omy_aromat280_mh") # removing triallelic snp


medLocus <- oxbow %>% calcHet() %>% arrange(expHet) %>% slice(ceiling(nrow(.)/2)) %>% pull(locus)
medLocusSNP <- oxbowSNP %>% calcHet() %>% arrange(expHet) %>% slice(ceiling(nrow(.)/2)) %>% pull(locus)

trem <- getLoci(oxbow)[getLoci(oxbow) != medLocus]
medData <- oxbow %>% removeLoci(trem)
oxbow_gma_mh <- exportGrandma(medData, baseline = TRUE)
gmaIn_mh <- createGmaInput(oxbow_gma_mh, perAlleleError = .005, dropoutProb = .005)

trem <- getLoci(oxbowSNP)[getLoci(oxbowSNP) != medLocusSNP]
medDataSNP <- oxbowSNP %>% removeLoci(trem)
oxbow_gma_snp <- exportGrandma(medDataSNP, baseline = TRUE)
gmaIn_snp <- createGmaInput(oxbow_gma_snp, perAlleleError = .005, dropoutProb = .005)

# edit missing data rates
gmaIn_mh$missingParams$Omy_RAD392622_mh.A1 <- c(3,97) # 97% genotyping success rate
gmaIn_snp$missingParams$Omy_ppie232_mh.A1 <- c(3,97) # 97% genotyping success rate

# edit allele frequencies to
# target exp het from Baetscher et al. 2018
gmaIn_mh$baselineParams$OmyOXBO20S$Omy_RAD392622_mh.A1 <- c(13,13,74)
gmaIn_snp$baselineParams$OmyOXBO20S$Omy_ppie232_mh.A1 <- c(12.5,87.5)

# now duplicate locus a particular number of times
# let's evaluate 100, 200, 300, and 400 loci

for(i in 1:399){
	# microhap
	n <- paste0(names(gmaIn_mh$genotypeErrorRates[1]), "_", i)
	gmaIn_mh$genotypeErrorRates <- c(gmaIn_mh$genotypeErrorRates, gmaIn_mh$genotypeErrorRates[1])
	names(gmaIn_mh$genotypeErrorRates)[i+1] <- n
	gmaIn_mh$genotypeKeys <- c(gmaIn_mh$genotypeKeys, gmaIn_mh$genotypeKeys[1])
	names(gmaIn_mh$genotypeKeys)[i+1] <- n
	gmaIn_mh$alleleKeys <- c(gmaIn_mh$alleleKeys, gmaIn_mh$alleleKeys[1])
	names(gmaIn_mh$alleleKeys)[i+1] <- n
	gmaIn_mh$missingParams <- c(gmaIn_mh$missingParams, gmaIn_mh$missingParams[1])
	names(gmaIn_mh$missingParams)[i+1] <- n

	gmaIn_mh$baselineParams$OmyOXBO20S <- c(gmaIn_mh$baselineParams$OmyOXBO20S, gmaIn_mh$baselineParams$OmyOXBO20S[1])
	names(gmaIn_mh$baselineParams$OmyOXBO20S)[i+1] <- n
	
	# snp
	n <- paste0(names(gmaIn_snp$genotypeErrorRates[1]), "_", i)
	gmaIn_snp$genotypeErrorRates <- c(gmaIn_snp$genotypeErrorRates, gmaIn_snp$genotypeErrorRates[1])
	names(gmaIn_snp$genotypeErrorRates)[i+1] <- n
	gmaIn_snp$genotypeKeys <- c(gmaIn_snp$genotypeKeys, gmaIn_snp$genotypeKeys[1])
	names(gmaIn_snp$genotypeKeys)[i+1] <- n
	gmaIn_snp$alleleKeys <- c(gmaIn_snp$alleleKeys, gmaIn_snp$alleleKeys[1])
	names(gmaIn_snp$alleleKeys)[i+1] <- n
	gmaIn_snp$missingParams <- c(gmaIn_snp$missingParams, gmaIn_snp$missingParams[1])
	names(gmaIn_snp$missingParams)[i+1] <- n
	
	gmaIn_snp$baselineParams$OmyOXBO20S <- c(gmaIn_snp$baselineParams$OmyOXBO20S, gmaIn_snp$baselineParams$OmyOXBO20S[1])
	names(gmaIn_snp$baselineParams$OmyOXBO20S)[i+1] <- n
	
	if(i == 99){
		gmaIn_mh_100 <- gmaIn_mh
		gmaIn_snp_100 <- gmaIn_snp
	}
	if(i == 199){
		gmaIn_mh_200 <- gmaIn_mh
		gmaIn_snp_200 <- gmaIn_snp
	}
	if(i == 299){
		gmaIn_mh_300 <- gmaIn_mh
		gmaIn_snp_300 <- gmaIn_snp
	}
}
gmaIn_mh_400 <- gmaIn_mh
gmaIn_snp_400 <- gmaIn_snp

print("false negative")

# false negative rates

fn_snp <- mclapply(list(gmaIn_snp_100, gmaIn_snp_200, gmaIn_snp_300, gmaIn_snp_400),
		 falseGrandma, relationship = "ssGP", llrToTest = seq(-10,20,1), seed = 7,
		 N = 10000, errorType = "falseNeg", mc.cores = numCores)
fn_mh <- mclapply(list(gmaIn_mh_100, gmaIn_mh_200, gmaIn_mh_300, gmaIn_mh_400),
		 falseGrandma, relationship = "ssGP", llrToTest = seq(-10,40,1), seed = 7,
		 N = 10000, errorType = "falseNeg", mc.cores = numCores)
fn_snp[[1]]
print("false positive strat")

# false positive unrel rates
system.time(
fp_unrel_snp <- mclapply(list(gmaIn_snp_100, gmaIn_snp_200, gmaIn_snp_300, gmaIn_snp_400),
					  falseGrandma, relationship = "ssGP", llrToTest = seq(-10,20,1), seed = 7,
					  itersPerMI = rep(1000000,20), errorType = "Unrel", mc.cores = numCores)
)
system.time(
fp_unrel_mh <- mclapply(list(gmaIn_mh_100, gmaIn_mh_200, gmaIn_mh_300, gmaIn_mh_400),
					 falseGrandma, relationship = "ssGP", llrToTest = seq(-10,40,1), seed = 7,
					 itersPerMI = rep(10000000,20), errorType = "Unrel", mc.cores = numCores)
)

print("false positive IS")

system.time(
	fp_unrel_snp_IS <- mclapply(list(gmaIn_snp_100, gmaIn_snp_200, gmaIn_snp_300, gmaIn_snp_400),
								  falseGrandma, relationship = "ssGP", llrToTest = seq(-10,20,1), seed = 7,
								  N = 1000000, errorType = "Unrel", method = "IS", mc.cores = numCores)
)
system.time(
	fp_unrel_mh_IS <- mclapply(list(gmaIn_mh_100, gmaIn_mh_200, gmaIn_mh_300, gmaIn_mh_400),
								 falseGrandma, relationship = "ssGP", llrToTest = seq(-10,40,1), seed = 7,
								 N = 1000000, errorType = "Unrel", method = "IS", mc.cores = numCores)
)

save.image("after_fn_unrel_lower.rda")
