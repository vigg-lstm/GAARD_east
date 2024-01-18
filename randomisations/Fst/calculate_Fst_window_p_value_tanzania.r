library(data.table)
library(stringr)
library(magrittr)

# Load the data
locations <- c('tanzania')
phenotype.column <- 'phenotype.moving.Fst'

all.fst.tables <- list()
for (l in locations){
	cat('Loading data for', l, '\n')
	load(paste('fst_randomisations_', l, '.Rdata', sep = ''))
	for (pop in names(fst.tables)){
		non.na.rows <- !is.na(fst.tables[[pop]][[phenotype.column]])
		all.fst.tables[[pop]] <- cbind(window.pos[[pop]][non.na.rows, ], fst.tables[[pop]][non.na.rows, ])
	}
}

fst.tables <- all.fst.tables
rm(all.fst.tables)

randomisation.columns <- grep('^r\\d{5}\\.', colnames(fst.tables[[1]]), value = T)
populations <- names(fst.tables)

# Calculate a P-value for each window based on the number of randomisations larger than the observed value
for (pop in populations){
	fst.tables[[pop]][, window.p := rowSums(apply(.SD, 2, function(x) x >= get(phenotype.column))) / ncol(.SD),
	                    .SDcols = randomisation.columns]
	fst.tables[[pop]][, log.p := -log10(window.p + 1/(length(..randomisation.columns)))]
}

# Identify outlier windows. Use the negative part of the distribution to tell us how far to look in the 
# positive part. 
first.pass.peaks <- function(fst.vector, modal.increments = 3, return.thresh = F){
	fst.modal.bin <- which.max(table(cut(fst.vector, breaks = seq(min(fst.vector, na.rm = T), max(fst.vector, na.rm = T), 0.0005))))
	fst.mode <- mean(as.numeric(strsplit(sub(']$', '', sub('^\\(', '', names(fst.modal.bin))), ',')[[1]]))
	neg.diff <- fst.mode - min(fst.vector, na.rm = T)
	thresh <- fst.mode + neg.diff * modal.increments
	if (return.thresh)
		thresh
	else
		fst.vector > thresh
}

for (pop in names(fst.tables))
	fst.tables[[pop]][, filter := first.pass.peaks(get(phenotype.column))]

p.threshold <- 0.01

# Next we want to use the randomised data to find which of those peaks we believe. 
# For each population, pull out the peaks. 
filtered.tables <- lapply(fst.tables, function(D) D[filter == T & window.p < p.threshold, ])

# Write a function to add transparency to a colour
add.transparency <- function(col, prop.alpha){
	rgb.val <- col2rgb(col)
	new.rgb <- rgb(rgb.val[1], rgb.val[2], rgb.val[3], max = 255, alpha = (1 - prop.alpha) * 255)
	new.rgb
}

chrom.sizes <- c('2R' = 61545105, '2L' = 49364325, '3R' = 53200684, '3L' = 41963435, 'X' = 24393108)

source('../../shared_functions/R_plotting.r')

# Show how these filters affect the peaks
plot.fst.filters <- function(fst.table, filename = NULL, num.randomisations = NULL, p.thresh = 0.01, plot.title = '', gaps = 5000000, filter.name = 'filter'){
	# Create the plot
	if (!missing(filename)){
		file.width = 6.5
		file.height = 3.5
		if (grepl('\\.eps', filename))
			postscript(filename, width = file.width, height = file.height, horizontal = F, onefile = FALSE, paper = "special")
		else if (grepl('\\.png', filename))
			png(filename, width = file.width, height = file.height, units = 'in', res = 600)
		else if (grepl('\\.tif', filename))
			tiff(filename, width = file.width, height = file.height, units = 'in', res = 600)
	}
	# Get vectors of start and end points for each chromosome (ie: cumulative sizes + gaps)
	ce <- cumsum(chrom.sizes + c(0, 0, gaps, 0, gaps))
	cs <- ce - chrom.sizes
	layout(matrix(c(rep(1,4),rep(2,1)), nrow = 5, ncol = 1))
	colours <- c(fst = add.transparency('orangered3', 0.2),
                 randomisations = add.transparency('grey50', 0.8))
	par(mar = c(0,4,1,2), mgp = c(2, 0.7, 0), family = 'Arial', xpd = NA) 
	num.randomisations <- ifelse(is.null(num.randomisations), 
	                             length(randomisation.columns),
	                             num.randomisations)
	# This odd way of getting a sequence is to make sure we get the right outcome if num.randomisations == 0
	r.columns <- randomisation.columns[seq(1, num.randomisations, length.out = num.randomisations)]
	fst.columns <- c(phenotype.column, r.columns)
	max.y <- max(c(max(fst.table[, ..fst.columns]), 0.05))
	min.y <- min(fst.table[, ..fst.columns])
	# Create the empty plot.
	plot(c(cs[1], ce[5]), c(min.y, max.y), xlim = c(cs[1] + gaps/2, ce[5] - gaps/2), type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', main = plot.title, cex.main = 1.2)
	
	# Get plotting positions for windows
	fst.table$genome.pos <- fst.table$window.pos + cs[as.character(fst.table$window.chrom)]
	# Add randomised data
	sapply(r.columns, function(x) fst.table[, lines(genome.pos, get(x), col = colours['randomisations'], lwd = 0.4), by = window.chrom])
	# Add true data
	fst.table[, lines(genome.pos, get(phenotype.column), col = colours['fst'], lwd = 1.2), by = window.chrom]
	# Add y axis
	fst.step.size <- ifelse(max.y > 0.2, 0.1, 0.05)
	axis(2, at = seq(0, max.y, fst.step.size))
	mtext('Fst', 2, 2, cex = 0.8)
	# Add peaks that pass filtering. Colour according to whether they have significant p-value
	if (length(p.thresh) == 1)
		p.thresh <- c(p.thresh, -1)
	p.colours <- c('orchid3', 'green', 'blue')
	filter.pass <- fst.table[[filter.name]]
	fst.table[filter.pass, points(genome.pos, get(phenotype.column), 
	                              pch = 21,
	                              bg = p.colours[(window.p < p.thresh[1]) + (window.p < p.thresh[2]) + 1], 
	                              col = colours['randomisations'], 
	                              cex = 1.1, lwd = .5)
	]
	
	# Now plot all chromosomes with, the position of each of the four detox gene regions and Ace1
	par(mar = c(1,4,0,2), mgp = c(2, 0.7, 0)) 
	add.chromosomes(gaps = gaps, chrom.sizes = chrom.sizes, gene.cex = 0.7, point.cex = 1, chrom.offset = -1.2, chrom.cex = 1.2)
	
	if (!missing(filename))
		dev.off()
}


for (pop in names(fst.tables))
	plot.fst.filters(fst.tables[[pop]], filename = paste(pop, 'peak_filter_plot.png', sep = '_'), plot.title = pop, p.thresh = p.threshold)


save.image('fst_filtered_windows_tanzania.Rdata')




