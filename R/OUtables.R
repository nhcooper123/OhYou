# Gavin Thomas April 2015
#
# Burrow to the working directory in your preferred way


	
# Set up empty table with 4 columns and 28 rows, assign names to columns
OUresTable <- as.data.frame(matrix(nrow=28, ncol=4))
colnames(OUresTable) <- c("Tree type", "Tree size", "Rejection rate", "Median alpha (95% HPD)")

OUresTable[1:7,"Tree type"] <- "b/d = 0"
OUresTable[8:14,"Tree type"] <- "b/d = 0.25"
OUresTable[15:21,"Tree type"] <- "b/d = 0.5"
OUresTable[22:28,"Tree type"] <- "b/d = 0.75"

OUresTable[seq(1, 28, 7), "Tree size"] <- 25
OUresTable[seq(2, 28, 7), "Tree size"] <- 50
OUresTable[seq(3, 28, 7), "Tree size"] <- 100
OUresTable[seq(4, 28, 7), "Tree size"] <- 150
OUresTable[seq(5, 28, 7), "Tree size"] <- 200
OUresTable[seq(6, 28, 7), "Tree size"] <- 500
OUresTable[seq(7, 28, 7), "Tree size"] <- 1000


# Now extract rejection rate and alpha estimates (median and 95% HPDs)

# Read in data with no error (rename data object to avoid confusion with same name of object in error sims)
load("Data/OU_simulations.rda")
noerror_res <- out

for (i in 1:28) {
	OUresTable[i,"Rejection rate"] <- sum(noerror_res[(1 + (i-1)*1000) : (i*1000), 5]>3.84)/1000	
	alpha <- quantile(noerror_res[(1 + (i-1)*1000) : (i*1000), 6], probs=c(0.5, 0.025, 0.975))	
	OUresTable[i,"Median alpha (95% HPD)"] <- paste(round(alpha[1], 3), " (", round(alpha[2], 3), "-", round(alpha[3], 3), ")", sep="")		
	}

write.csv(OUresTable, "Manuscript/tables/OU_bd_noerror.csv")




# Read in data with error (rename data object to avoid confusion with same name of object in error sims)
load("Data/OU_error_simulations.rda")
error_res <- out

OUerrorTable <- as.data.frame(matrix(nrow=7, ncol=8))
colnames(OUerrorTable) <- c("Rejection (0% error)", "Rejection (1% error)", "Rejection (5% error)", "Rejection (10% error)", "Alpha (0% error)", "Alpha (1% error)", "Alpha (5% error)", "Alpha (10% error)")

rownames(OUerrorTable)<- c(25,50,100,150,200,500,1000)

for (i in 1:7) {
	OUerrorTable[i,"Rejection (0% error)"] <- sum(noerror_res[(1 + (i-1)*1000) : (i*1000), 5]>3.84)/1000	
	alpha <- quantile(noerror_res[(1 + (i-1)*1000) : (i*1000), 6], probs=c(0.5, 0.025, 0.975))	
	OUerrorTable[i,"Alpha (0% error)"] <- paste(round(alpha[1], 3), " (", round(alpha[2], 3), "-", round(alpha[3], 3), ")", sep="")		
	OUerrorTable[i,2] <- sum(error_res[(1 + (i-1)*1000) : (i*1000), 5]>3.84)/1000
	alpha <- quantile(error_res[(1 + (i-1)*1000) : (i*1000), 6], probs=c(0.5, 0.025, 0.975))	
	OUerrorTable[i,6] <- paste(round(alpha[1], 3), " (", round(alpha[2], 3), "-", round(alpha[3], 3), ")", sep="")

	OUerrorTable[i,3] <- sum(error_res[((1 + (i-1)*1000)+7000) : ((i*1000)+7000), 5]>3.84)/1000
	alpha <- quantile(error_res[((1 + (i-1)*1000)+7000) : ((i*1000)+7000), 6], probs=c(0.5, 0.025, 0.975))	
	OUerrorTable[i,7] <- paste(round(alpha[1], 3), " (", round(alpha[2], 3), "-", round(alpha[3], 3), ")", sep="")

	OUerrorTable[i,4] <- sum(error_res[((1 + (i-1)*1000)+14000) : ((i*1000)+14000), 5]>3.84)/1000
	alpha <- quantile(error_res[((1 + (i-1)*1000)+14000) : ((i*1000)+14000), 6], probs=c(0.5, 0.025, 0.975))	
	OUerrorTable[i,8] <- paste(round(alpha[1], 3), " (", round(alpha[2], 3), "-", round(alpha[3], 3), ")", sep="")
		}

write.csv(OUerrorTable, "Manucript/tables/OU_yule_error.csv")


# Now for rate variable through time trees
# Set up empty table with 4 columns and 28 rows, assign names to columns
OUresVarRateTable <- as.data.frame(matrix(nrow=28, ncol=4))
colnames(OUresVarRateTable) <- c("Tree type", "Tree size", "Rejection rate", "Median alpha (95% HPD)")

OUresVarRateTable[1:7,"Tree type"] <- "Slow speed-up"
OUresVarRateTable[8:14,"Tree type"] <- "Rapid speed-up"
OUresVarRateTable[15:21,"Tree type"] <- "Slow slow-down"
OUresVarRateTable[22:28,"Tree type"] <- "Rapid slow-down"

OUresVarRateTable[seq(1, 28, 7), "Tree size"] <- 25
OUresVarRateTable[seq(2, 28, 7), "Tree size"] <- 50
OUresVarRateTable[seq(3, 28, 7), "Tree size"] <- 100
OUresVarRateTable[seq(4, 28, 7), "Tree size"] <- 150
OUresVarRateTable[seq(5, 28, 7), "Tree size"] <- 200
OUresVarRateTable[seq(6, 28, 7), "Tree size"] <- 500
OUresVarRateTable[seq(7, 28, 7), "Tree size"] <- 1000


# Now extract rejection rate and alpha estimates (median and 95% HPDs)

for (i in 1:28) {
	OUresVarRateTable[i,"Rejection rate"] <- sum(noerror_res[((1 + (i-1)*1000)+28000) : ((i*1000)+28000), 5]>3.84)/1000	
	alpha <- quantile(noerror_res[((1 + (i-1)*1000)+28000) : ((i*1000)+28000), 6], probs=c(0.5, 0.025, 0.975))	
	OUresVarRateTable[i,"Median alpha (95% HPD)"] <- paste(round(alpha[1], 3), " (", round(alpha[2], 3), "-", round(alpha[3], 3), ")", sep="")		
	}

write.csv(OUresVarRateTable, "Manuscript/tables/OU_var_rate_noerror.csv")



####################
## Bayestraits tables
OUBayes <- read.csv("Data/MCMC_ExpLam10.csv", stringsAsFactors=FALSE)


OUBayesTable <- as.data.frame(matrix(nrow=28, ncol=4))
colnames(OUBayesTable) <- c("Tree type", "Tree size", "Rejection rate", "Median alpha (95% HPD)")

OUBayesTable[1:7,"Tree type"] <- "b/d = 0"
OUBayesTable[8:14,"Tree type"] <- "b/d = 0.25"
OUBayesTable[15:21,"Tree type"] <- "b/d = 0.5"
OUBayesTable[22:28,"Tree type"] <- "b/d = 0.75"

OUBayesTable[seq(1, 28, 7), "Tree size"] <- 25
OUBayesTable[seq(2, 28, 7), "Tree size"] <- 50
OUBayesTable[seq(3, 28, 7), "Tree size"] <- 100
OUBayesTable[seq(4, 28, 7), "Tree size"] <- 150
OUBayesTable[seq(5, 28, 7), "Tree size"] <- 200
OUBayesTable[seq(6, 28, 7), "Tree size"] <- 500
OUBayesTable[seq(7, 28, 7), "Tree size"] <- 1000


nms <- c("yule25", "yule50", "yule100", "yule150", "yule200", "yule500", "yule1000", "bdLow25", "bdLow50", "bdLow100", "bdLow150", "bdLow200", "bdLow500", "bdLow1000", "bdMid25", "bdMid50", "bdMid100", "bdMid150", "bdMid200", "bdMid500", "bdMid1000", "bdHigh25", "bdHigh50", "bdHigh100", "bdHigh150", "bdHigh200", "bdHigh500", "bdHigh1000")

for (i in 1:28) {
x <- OUBayes[OUBayes[,"Model"] == paste(nms[i]),]
OUBayesTable[i,"Rejection rate"] <- round(sum(x[,"Exp_L10_2Bfactor"] > 2, na.rm=T)/sum(!is.na(x[,"Exp_L10_2Bfactor"])), 3)

alpha <- quantile(x[,"Exp_L10_Alpha_Mode"], probs=c(0.5, 0.025, 0.975))

    OUBayesTable[i,"Median alpha (95% HPD)"] <- paste(round(alpha[1], 3), " (", round(alpha[2], 3), "-", round(alpha[3], 3), ")", sep="")

}


write.csv(OUBayesTable, "Manuscript/tables/OU_bd_noerror_BayesTraits.csv")



## BayesTraits with error - data are in the same table as no error simulations


OUBayeserrorTable <- as.data.frame(matrix(nrow=7, ncol=8))
colnames(OUBayeserrorTable) <- c("Rejection (0% error)", "Rejection (1% error)", "Rejection (5% error)", "Rejection (10% error)", "Alpha (0% error)", "Alpha (1% error)", "Alpha (5% error)", "Alpha (10% error)")

rownames(OUBayeserrorTable)<- c(25,50,100,150,200,500,1000)

nms <- c("yule25", "yule50", "yule100", "yule150", "yule200", "yule500", "yule1000")

for (i in 1:7) {
    
    x1 <- OUBayes[OUBayes[,"Model"] == paste(nms[i]),]
    x2 <- OUBayes[OUBayes[,"Model"] == paste(nms[i], "_err1", sep=""),]
    x3 <- OUBayes[OUBayes[,"Model"] == paste(nms[i], "_err5", sep=""),]
    x4 <- OUBayes[OUBayes[,"Model"] == paste(nms[i], "_err10", sep=""),]
    
    OUBayeserrorTable[i,"Rejection (0% error)"] <- round(sum(x1[,"Exp_L10_2Bfactor"] > 2, na.rm=T)/sum(!is.na(x[,"Exp_L10_2Bfactor"])), 3)
    
    OUBayeserrorTable[i,"Rejection (1% error)"] <- round(sum(x2[,"Exp_L10_2Bfactor"] > 2, na.rm=T)/sum(!is.na(x[,"Exp_L10_2Bfactor"])), 3)
        
    OUBayeserrorTable[i,"Rejection (5% error)"] <- round(sum(x3[,"Exp_L10_2Bfactor"] > 2, na.rm=T)/sum(!is.na(x[,"Exp_L10_2Bfactor"])), 3)
    
    OUBayeserrorTable[i,"Rejection (10% error)"] <- round(sum(x4[,"Exp_L10_2Bfactor"] > 2, na.rm=T)/sum(!is.na(x[,"Exp_L10_2Bfactor"])), 3)

    alpha <- quantile(x1[,"Exp_L10_Alpha_Mode"], probs=c(0.5, 0.025, 0.975))
    OUBayeserrorTable[i,"Alpha (0% error)"] <- paste(round(alpha[1], 3), " (", round(alpha[2], 3), "-", round(alpha[3], 3), ")", sep="")
    
    alpha <- quantile(x2[,"Exp_L10_Alpha_Mode"], probs=c(0.5, 0.025, 0.975))
    OUBayeserrorTable[i,"Alpha (1% error)"] <- paste(round(alpha[1], 3), " (", round(alpha[2], 3), "-", round(alpha[3], 3), ")", sep="")
    
    alpha <- quantile(x3[,"Exp_L10_Alpha_Mode"], probs=c(0.5, 0.025, 0.975))
    OUBayeserrorTable[i,"Alpha (5% error)"] <- paste(round(alpha[1], 3), " (", round(alpha[2], 3), "-", round(alpha[3], 3), ")", sep="")
    
    alpha <- quantile(x4[,"Exp_L10_Alpha_Mode"], probs=c(0.5, 0.025, 0.975))
    OUBayeserrorTable[i,"Alpha (10% error)"] <- paste(round(alpha[1], 3), " (", round(alpha[2], 3), "-", round(alpha[3], 3), ")", sep="")
}

write.csv(OUBayeserrorTable, "Manuscript/tables/OU_yule_error_BayesTraits.csv")




# Now for rate variable through time trees for BayesTraits output
# Set up empty table with 4 columns and 28 rows, assign names to columns
OUBayesresVarRateTable <- as.data.frame(matrix(nrow=28, ncol=4))
colnames(OUBayesresVarRateTable) <- c("Tree type", "Tree size", "Rejection rate", "Median alpha (95% HPD)")

OUBayesresVarRateTable[1:7,"Tree type"] <- "Slow speed-up"
OUBayesresVarRateTable[8:14,"Tree type"] <- "Rapid speed-up"
OUBayesresVarRateTable[15:21,"Tree type"] <- "Slow slow-down"
OUBayesresVarRateTable[22:28,"Tree type"] <- "Rapid slow-down"

OUBayesresVarRateTable[seq(1, 28, 7), "Tree size"] <- 25
OUBayesresVarRateTable[seq(2, 28, 7), "Tree size"] <- 50
OUBayesresVarRateTable[seq(3, 28, 7), "Tree size"] <- 100
OUBayesresVarRateTable[seq(4, 28, 7), "Tree size"] <- 150
OUBayesresVarRateTable[seq(5, 28, 7), "Tree size"] <- 200
OUBayesresVarRateTable[seq(6, 28, 7), "Tree size"] <- 500
OUBayesresVarRateTable[seq(7, 28, 7), "Tree size"] <- 1000


nms <- c("slowIncrease25", "slowIncrease50", "slowIncrease100", "slowIncrease150", "slowIncrease200", "slowIncrease500", "slowIncrease1000", "rapidIncrease25", "rapidIncrease50", "rapidIncrease100", "rapidIncrease150", "rapidIncrease200", "rapidIncrease500", "rapidIncrease1000", "rapidDecrease25", "rapidDecrease50", "rapidDecrease100", "rapidDecrease150", "rapidDecrease200", "rapidDecrease500", "rapidDecrease1000", "rapidIncrease25", "rapidIncrease50", "rapidIncrease100", "rapidIncrease150", "rapidIncrease200", "rapidIncrease500", "rapidIncrease1000")

# Now extract rejection rate and alpha estimates (median and 95% HPDs)

for (i in 1:28) {
    
    x <- OUBayes[OUBayes[,"Model"] == paste(nms[i]),]
    OUBayesresVarRateTable[i,"Rejection rate"] <- round(sum(x[,"Exp_L10_2Bfactor"] > 2, na.rm=T)/sum(!is.na(x[,"Exp_L10_2Bfactor"])), 3)
    
    alpha <- quantile(x[,"Exp_L10_Alpha_Mode"], probs=c(0.5, 0.025, 0.975))
    
    OUBayesresVarRateTable[i,"Median alpha (95% HPD)"] <- paste(round(alpha[1], 3), " (", round(alpha[2], 3), "-", round(alpha[3], 3), ")", sep="")
}

write.csv(OUBayesresVarRateTable, "Manuscript/tables/OU_var_rate_noerror_BayesTraits.csv")


