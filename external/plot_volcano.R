# PLOT VOLCANO # by Vincent

# This function takes a data frame with three columns and plots a volcano plot!

# Submit data in the form of:

# Column 1: data label text
# Column 2: Log fold changes
# Column 3: p values
# Column 4: base means

# OR submit a DESeq results file and specify the argument: input.format="DESeq"

# Note: colour input has to be in the form of RGB vectors. i.e: c(255, 0, 0) for pure red.

plot.volcano <- function(data, input.format, colour1, colour2, fold.change, p.level, basemean.cut, text, fade_style, lines, plot.legend, legendPosition, cex,...) {

# Input:
# data = dataframe(dataLabels,logFoldChange,Pvalue,BaseMean) (or see below to use a DESeq results object).

# input.format = can be specified to be a "DESeq" format. For ease of use...

# colour1 = rgb(r1, g1, b1, o1)

# colour2 = rgb(r2, g2, b2, o2)

# fold.change = FC cut off, (non-log) This is log2d in the function.

# p.level = p value cutoff  (non-log)

# basemean.cut = base mean cut off. This is the 1st quartile by default.

# text = (T / F) display data labels under data points

# fade_style = (1 / 2) 1 is fold-change only (default). 2 includes cross-fade for p.values.

# lines = (T / F) includes the cut off lines as dashed lines.

if(missing(input.format)){input.format = "df"}
if(missing(colour1)){colour1 = c(200, 80, 80)}
if(missing(colour2)){colour2 = c(80, 200, 80)}
if(missing(fold.change)){fold.change = 2}
if(missing(p.level)){p.level = 0.05}
if(missing(basemean.cut)){basemean.cut = summary(data[,4])[2]}
if(missing(text)){do_you_want_text = F} else {do_you_want_text = T}
if(missing(fade_style)){fade_style = 1}
if(missing(lines)){do_you_want_lines = F} else {do_you_want_lines = T}
if(missing(plot.legend)){plot.legend = T}
if(missing(legendPosition)){legendPosition = "topleft"}
if(missing(cex)){cex=2}
if(input.format == "DESeq") {data = data.frame(rownames(data), data[,2], data[,6], data[,1])}

# Data manipulation
	
	data1 <- data[data[,2] > 0,]
	# check for positive fold change genes - create dummy otherwise
	if(length(data1[,1]) != 0) {pos = data1} else {pos <- data.frame(0,0,1,0)}
	pos <- pos[which(as.character(pos[,3]) != "NA"),]

	data2 <- data[data[,2] < 0,]
	# check for negative fold change genes - create dummy otherwise
	if(length(data2[,1]) != 0) {neg = data2} else {neg <- data.frame(0,0,1,0)}
	neg <- neg[which(as.character(neg[,3]) != "NA"),]

	
# remove zeros
	pos[pos[,3] == 0,3] <- min(pos[pos[,3] != 0,3])
	neg[neg[,3] == 0,3] <- min(neg[neg[,3] != 0,3])

	pos[,3] <- -log10(pos[,3])
	neg[,3] <- -log10(neg[,3])

# Let's go!
#####
#	pos <- data1
#	neg <- data2
	#if(length(data1[,1]) != 0) {pos <- data1} else {pos <- data.frame(0,0,0,0)}
	#if(length(data2[,1]) != 0) {neg <- data2} else {neg <- data.frame(0,0,0,0)}

	n.col = colour1
	p.col = colour2

# Initialise plot
	par(cex=cex); plot(c(min(neg[,2])+min(neg[,2])/15,max(pos[,2])+max(pos[,2])/15),c(0,max(neg[,3],pos[,3])+max(neg[,3],pos[,3])/15), pch=NA, col="white", xlab="log2 fold-change", ylab="-log10 adjusted p-value",...)
#par(cex=cex); plot(c(min(neg[,2])+min(neg[,2])/15,max(pos[,2])+max(pos[,2])/15),c(0,max(neg[,3],pos[,3])+max(neg[,3],pos[,3])/15), pch=NA, col="white", if(!exists("xlab")){xlab="log 2 fold-change"}, if(!exists("ylab")){ylab="-log 10 P"}, ...)

if(length(data2[,1]) != 0) {
# negative change
	x <- neg[,2]
	y <- neg[,3]
	z <- neg[,4]
	res.vec <- which(y > -log10(p.level) & x < -log2(fold.change) & z > basemean.cut)

# text?
	if(do_you_want_text == T){if(length(x[res.vec]) != 0){text(x[res.vec], y[res.vec], neg[res.vec,1])}}

# colour grading
			x_r2 <- x/(min(x)/2)
			x_r2[x_r2 > 0.6] <- 0.6
			x_r <- x/min(x) 
			y_r <- y/max(y)
			z_r <- x_r * y_r
			z_r <- z_r/max(z_r)
			if(fade_style == 1){r = x_r} else {r = z_r}
			red <- 1-(r*(1-n.col[1]/255))
			green <- 1-(r*(1-n.col[2]/255))
			blue <- 1-(r*(1-n.col[3]/255))
			cols <- rgb(red,green,blue)
# plot all
			points(x, y, pch=16, col=cols)

# plot significant only
			points(x[res.vec], y[res.vec], pch=21, bg=rgb(n.col[1]/255,n.col[2]/255,n.col[3]/255,1), col="black", lwd=2)
}

if(length(data1[,1]) != 0) {
# positive change
	x <- pos[,2]
	y <- pos[,3]
	z <- pos[,4]
	res.vec <- which(y > -log10(p.level) & x > log2(fold.change) & z > basemean.cut)

# text?
	if(do_you_want_text == T){if(length(x[res.vec]) != 0) {text(x[res.vec], y[res.vec], pos[res.vec,1])}}

# colour grading
			x_r2 <-  x/(max(x)/2)
			x_r2[x_r2 > 0.6] <- 0.6
			x_r <- x/max(x)
			y_r <- y/max(y)
			z_r <- x_r * y_r
			z_r <- z_r/max(z_r)
			if(fade_style == 1){r = x_r} else {r = z_r}
			red <- 1-(r*(1-p.col[1]/255))
			green <- 1-(r*(1-p.col[2]/255))
			blue <- 1-(r*(1-p.col[3]/255))
			cols <- rgb(red,green,blue)
# plot all
			points(x, y, pch=16, col=cols)

# plot significant only
			points(x[res.vec], y[res.vec], pch=21, bg=rgb(p.col[1]/255,p.col[2]/255,p.col[3]/255,1), col="black", lwd=2)

# lines?
	if(do_you_want_lines == T){
	lines(c(-10000,10000),rep(-log10(p.level),2), lty=2, lwd=1)
	lines(c(-log2(fold.change),-log2(fold.change)),c(-10000,10000), lty=2, lwd=1)
	lines(c(log2(fold.change),log2(fold.change)),c(-10000,10000), lty=2, lwd=1)
	}

# Legend?

	if(plot.legend == T){
	par(lend=2); legend(legendPosition, legend=c("Increased in CF", "Increased in EPI"), pt.bg=c(rgb(n.col[1]/255,n.col[2]/255,n.col[3]/255,1), rgb(p.col[1]/255,p.col[2]/255,p.col[3]/255,1)), lty=0, pch=21, lwd=1, bty="n")
	}
}

##################
# function(data, colour1, colour2, fold.change, p.level, do_you_want_text, fade_style, do_you_want_lines, plot.legend, legendPosition) 
# Input:
# data = dataframe(dataLabels,logFoldChange,Pvalue,BaseMean) (or see below to use a DESeq results object).

# input.format = can be specified to be a "DESeq" format. For ease of use...

# colour1 = rgb(r1, g1, b1, o1)

# colour2 = rgb(r2, g2, b2, o2)

# fold.change = FC cut off, (non-log) This is log2d in the function.

# p.level = p value cutoff  (non-log)

# basemean.cut = base mean cut off. This is the 1st quartile by default.

# text = (T / F) display data labels under data points

# fade_style = (1 / 2) 1 is fold-change only (default). 2 includes cross-fade for p.values.

# lines = (T / F) includes the cut off lines as dashed lines.

}

############## EXAMPLE ##############################
# Some dodgy distributions here...
	paste(LETTERS[round(runif(5, 1,26))],collapse="")
		data.labels <- sapply(1:1000, function(x) {paste(LETTERS[round(runif(5, 1,26))],collapse="")})
		logfoldchange <- rnorm(1000,0,2)
		pvalues <- abs(rnorm(1000,0.3,0.5))^2
			pvalues[pvalues > 1] <- 1
		basemeans <- abs(rnorm(1000,500,300))

# Create dataframe with the correct column order...
	t1 <- data.frame(data.labels, 
			logfoldchange,
			pvalues,
			basemeans
			)

# Plot multiple versions...
	par(mfrow=c(2,2), mar=c(5,5,2,2));

# Here is the function:
	plot.volcano(t1, cex=1.5)
	plot.volcano(t1, legend=F, cex=1.5)
	plot.volcano(t1, legend=F, lines=T, cex=1.5, colour1=c(195,65,173), colour2=c(65,220,205), fade=1)
	plot.volcano(t1, legend=F, lines=T, text=T, cex=1.5, fade=2); par(mfrow=c(1,1))
######################################################

