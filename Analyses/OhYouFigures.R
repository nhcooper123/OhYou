# Making figures 

# Read in data
ds <- read.delim("Data/papers.txt")

# Calculate %
ds$eco <- (ds$ecology / ds$ecology.total) * 100
ds$evo <- (ds$evolution / ds$evolution.total) * 100
ds$pal <- (ds$paleo / ds$paleo.total) * 100

# Plot figure 1
pdf("Manuscript/Figures/OhYou_Figure1.pdf")
par(bty = "l")
plot(evo ~ Year, data = ds, las = 1, pch = 16, cex = 1.2, cex.lab = 1.2, cex.axis = 1.2,
     ylab = "OU publications (% of total papers)", xlab = "Year published", ylim = c(0, 0.7))
points(eco ~ Year, data = ds, cex = 1.2)
points(pal ~ Year, data = ds, cex = 1.2, pch = 6)
legend("topleft", legend = c("ecology", "evolution", "palaeontology"), pch = c(1, 16, 6), cex = 1.2, bty = "n")
dev.off()

# Read in data for table 2 and figure 3 and reshape
ds1 <- read.delim("Data/literature.txt")
ds1 <- subset(ds1, Noptima != "?")
ds1$Noptima <- ds1$Noptima[, drop = TRUE]

# Table 2
omit <- which(ds1$General.use == "other" | ds1$General.use == "convergent evolution" |
	  ds1$General.use == "ancestral state reconstruction" | 
	  ds1$General.use == "phylogenetic inertia" | ds1$General.use == "phylogenetic signal" )

ds2 <- ds1[-omit, ]
ds2$General.use <- ds2$General.use[, drop = TRUE]

summary.table <- tapply(ds2$Paper, list(ds2$General.use, ds2$Noptima), length)

# Figure 3
pdf("Manuscript/Figures/OhYou_Figure3.pdf")
  par(bty = "l")
  hist(ds1$Ntips[which(ds1$Ntips <= 1000)], las = 1, cex.lab = 1.2, cex.axis = 1.2,
       ylab = "Number of papers", xlab = "Number of taxa", ylim = c(0, 30), main = "",
       breaks = seq(from = 0, to = 800, by = 25))
dev.off()
