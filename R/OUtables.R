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
load("OU_simulations.rda")
noerror_res <- out

for (i in 1:28) {
	OUresTable[i,"Rejection rate"] <- sum(noerror_res[(1 + (i-1)*1000) : (i*1000), 5]>3.84)/1000	
	alpha <- quantile(noerror_res[(1 + (i-1)*1000) : (i*1000), 6], probs=c(0.5, 0.025, 0.975))	
	OUresTable[i,"Median alpha (95% HPD)"] <- paste(round(alpha[1], 3), " (", round(alpha[2], 3), "-", round(alpha[3], 3), ")", sep="")		
	}

write.csv(OUresTable, "tables/OU_bdTrees.csv")




# Read in data with error (rename data object to avoid confusion with same name of object in error sims)
load("OU_error_simulations.rda")
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


write.csv(OUerrorTable, "tables/OU_yule_error_Trees.csv")



















# Set up empty table with 10 columns and 28 rows, assign names to columns
OUresTable <- as.data.frame(matrix(nrow=28, ncol=10))
colnames(OUresTable) <- c("Tree type", "Tree size", "Rejection (0% error)", "Rejection (1% error)", "Rejection (5% error)", "Rejection (10% error)", "Alpha (0% error)", "Alpha (1% error)", "Alpha (5% error)", "Alpha (10% error)")

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
load("OU_simulations.rda")
noerror_res <- out

for (i in 1:28) {
	OUresTable[i,"Rejection (0% error)"] <- sum(noerror_res[(1 + (i-1)*1000) : (i*1000), 5]>3.84)/1000	
	alpha <- quantile(noerror_res[(1 + (i-1)*1000) : (i*1000), 6], probs=c(0.5, 0.025, 0.975))	
	OUresTable[i,"Alpha (0% error)"] <- paste(round(alpha[1], 3), " (", round(alpha[2], 3), "-", round(alpha[3], 3), ")", sep="")		
	}


# Read in data with error (rename data object to avoid confusion with same name of object in error sims)
load("OU_error_simulations.rda")
error_res <- out

for (i in 1:28) {
	OUresTable[i,"Rejection (0% error)"] <- sum(noerror_res[(1 + (i-1)*1000) : (i*1000), 5]>3.84)/1000	
	alpha <- quantile(noerror_res[(1 + (i-1)*1000) : (i*1000), 6], probs=c(0.5, 0.025, 0.975))	
	OUresTable[i,"Alpha (0% error)"] <- paste(round(alpha[1], 3), " (", round(alpha[2], 3), "-", round(alpha[3], 3), ")", sep="")		
	}

