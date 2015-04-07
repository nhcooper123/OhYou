# Making figures 

# Read in data
ds <- read.delim("Data/papers.txt")

# Calculate %
ds$eco <- (ds$ecology / ds$ecology.total) * 100
ds$evo <- (ds$evolution / ds$evolution.total) * 100
ds$pal <- (ds$paleo / ds$paleo.total) * 100

# Plot figure 1

png("Manuscript/Figures/PapersThruTime.png")
par(bty = "l")
plot(evo ~ Year, data = ds, las = 1, pch = 16, cex = 1.2, cex.lab = 1.2, cex.axis = 1.2,
     ylab = "OU publications (% of total papers)", xlab = "Year published", ylim = c(0, 0.7))
points(eco ~ Year, data = ds, cex = 1.2)
points(pal ~ Year, data = ds, cex = 1.2, pch = 6)
legend("topleft", legend = c("ecology", "evolution", "palaeontology"), pch = c(1, 16, 6), cex = 1.2, bty = "n")
dev.off()