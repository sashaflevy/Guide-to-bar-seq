
setwd("~/OneDrive - Leland Stanford Junior University/Github/Guide-to-bar-seq/")

#Functions
randDNA = function(n){
 sample(c("A","C","T","G"), n, replace=TRUE)
}

hamming = function(test, query = test){
  sum(query != test)
}

nearest.neighbors = function(bc.length, bc.number){
  #Make data matrix
  bc = matrix(NA, bc.number, bc.length)
  for(i in 1:bc.number){
    bc[i,] = randDNA(bc.length)
  }
  
  #Find Nearest Neighbor Distance
  d = 1:bc.number
  for(i in 1:bc.number){
    d[i] = min(apply(bc, 1, hamming, query = bc[i,])[-i])
  }
  return(d)
}

#Import data
nnfile  = "nearest_neighbor100length10_30.csv"
x2 = as.matrix(read.csv(nnfile, header = F))
rownames(x2) = 10:30
nnfile  = "nearest_neighbor1000length10_30.csv"
x3 = as.matrix(read.csv(nnfile, header = F))
rownames(x3) = 10:30
nnfile  = "nearest_neighbor10000length10_30.csv"
x4 = as.matrix(read.csv(nnfile, header = F))
rownames(x4) = 10:30
nnfile  = "nearest_neighbor100000length10_30.csv"
x5 = as.matrix(read.csv(nnfile, header = F))
rownames(x5) = 10:30

min2 = apply(x2, 1, min)
min3 = apply(x3, 1, min)
min4 = apply(x4, 1, min)
min5 = apply(x5, 1, min)
plot(names(min4), min4, ylab = "Nearest barcodes (Hammming distance)", xlab = "Barcode length", pch = 19, col = rgb(0.7,0,1,0.3), ylim = c(0,15), cex.lab = 1.5)
abline(lm(min4 ~ as.numeric(names(min4))), lwd = 2, col = rgb(0.7,0,1))
points(names(min3), min3, ylab = "Hamming distance betwwen the two nearest barcodes", xlab = "Barcode length", pch = 19, col = rgb(0.4,.4,0.7,0.3), ylim = c(0,15))
abline(lm(min3 ~ as.numeric(names(min3))), lwd = 2, col = rgb(0.4,.4,0.7))
points(names(min5), min5, ylab = "Hamming distance betwwen the two nearest barcodes", xlab = "Barcode length", pch = 19, col = rgb(0,0.5,1,0.3), ylim = c(0,15))
abline(lm(min5 ~ as.numeric(names(min5))), lwd = 2, col = rgb(0,0.5,1))
abline(h = 4, col = "grey", lwd = 4, lty = 5)
text(28, 4.5, "4 mismatches", col = "grey", cex = 1.2)

legend(10, 15, legend = c(as.expression(bquote(10^3)), as.expression(bquote(10^4)), as.expression(bquote(10^5))), lwd = 2, col=c(  rgb(0.4,.4,0.7), rgb(0.7,0,1) , rgb(0,0.5,1)), title = "Library size")



