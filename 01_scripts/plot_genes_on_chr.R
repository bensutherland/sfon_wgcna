# Plot genes on chromosomes
# input: chromosome lengths, bed file, transcript names of interest
par(mfrow=c(1,1), mar = c(4,4,3,3), cex = 0.6)

# use plot to draw empty graph
plot(1, type="n", xlab="", ylab="", xlim=c(1, 10), ylim=c(0, 10), yaxt = "n")
## note: replace xlim with the max length of chromosomes of interest
axis(2, at=1:10, labels=paste("chr_", seq(1:10)), las = 1)


# Plot multiple chromosomes this way:
segments(x0 = 1, y0 = 10, x1 = 10, y1 = 10, col = "black", lwd = 1)

# add genes this way
chr.pos <- 10 # this would be the level that the chr is at on the graph
start <- 2 # this would be the starting position of the gene
stop <- 3 # this would be the ending position of the gene

segments(x0 = start, y0 = chr.pos, x1 = stop, y1 = chr.pos, col = "blue", lwd = 5)


# Plot a couple examples..