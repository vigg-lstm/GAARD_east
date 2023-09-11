library(data.table)
library(magrittr)
library(ape)
library(stringr)

h12.tables <- readRDS('h12_filtered_windows_tanzania.RDS')
study.pops <- names(h12.tables)
# For some reason that I don't understand, these tables appear to have been loaded as immutable, such that 
# if you try to add a column (eg: h12.tables[[1]]$diff[, test := 1] ), it doesn't work. I have no idea why
# this is. To fix it, we basically need to recreate the object by copying the data from its original form. 
# Very odd.
for (pop in study.pops){
	for (table.type in names(h12.tables[[pop]])){
				h12.tables[[pop]][[table.type]] <- copy(h12.tables[[pop]][[table.type]])
	}
}

# For each table, identify regions around peaks, which we define as a window that is a peak, or within 50 
# windows of one.
peak.regions <- function(x, region.buffer = 20, value.to.check = T) {
	regions <- sapply(1:length(x), 
	                  function(i) 
	                  	any(x[max(1, i-region.buffer):min(length(x), i+region.buffer)] == value.to.check))
	region.starts <- (regions == T) & (c(F, regions[-length(regions)]) == F)
	region.codes <- cumsum(region.starts)
	region.codes[regions == F] <- 0
	region.codes
}

for (pop in study.pops){
	h12.tables[[pop]]$diff[, peak.region := peak.regions(is.peak)]
}

# Now pull out the peak regions separately
randomisation.ids <- names(h12.tables[[1]]$diff) %>%
                     grep('r\\d{5}', ., value = T)
p.threshold = 0.01

# Function to plot the H12 data
source('../../shared_functions/R_plotting.r')
chrom.sizes <- c('2R' = 61545105, '2L' = 49364325, '3R' = 53200684, '3L' = 41963435, 'X' = 24393108)
gff.path <- '../../data/VectorBase-57_AgambiaePEST.gff'
gff <- as.data.table(read.gff(gff.path, GFF3 = T))
gff[, seqid := sub('AgamP4_', '', seqid)]

plot.h12.diff.region <- function(h12.table, 
                                 filename = NULL, 
                                 num.randomisations = NULL, 
                                 chrom = NULL, 
                                 p.thresh = 0.01, 
                                 plot.title = '', 
                                 filter.name = 'is.peak'){
	if (is.null(chrom)){
		chrom <- unique(h12.table$chromosome)
		if (length(chrom) != 1)
			stop('Region should include exactly one chromosome')
	}
		
	start.pos <- min(h12.table$startpoint)
	end.pos <- max(h12.table$endpoint)
	
	# Get the list of genes present in the window
	gene.gff <- gff[seqid == chrom & 
	                start < end.pos & 
	                end > start.pos &
	                type == 'protein_coding_gene', ]
	gene.gff[, gene.id := unique(str_extract(attributes, 'AGAP\\d{6}'))]
	gene.gff[, gene.name := str_extract(attributes, '(?<=Name=)[^;]+') %>%
	                        str_to_title]
	gene.gff$gene.name[is.na(gene.gff$gene.name)] <- gene.gff$gene.id[is.na(gene.gff$gene.name)]
	# Now the list of Exons
	exon.gff <- gff[seqid == chrom & 
	                start < end.pos & 
	                end > start.pos &
	                type == 'exon', ]
	
	# Create the plot
	if (!is.null(filename)){
		file.width = 6.5
		file.height = 3.5
		if (grepl('\\.eps', filename))
			postscript(filename, width = file.width, height = file.height, horizontal = F, onefile = FALSE, paper = "special")
		else if (grepl('\\.png', filename))
			png(filename, width = file.width, height = file.height, units = 'in', res = 600)
		else if (grepl('\\.tif', filename))
			tiff(filename, width = file.width, height = file.height, units = 'in', res = 600)
	}
	layout(matrix(c(rep(1,4),rep(2,1)), nrow = 5, ncol = 1))
	colours <- c(h12 = lighten.col('limegreen', alpha = 0.8),
                 randomisations = lighten.col('grey50', alpha = 0.15))
	par(mar = c(0,4,1,1), mgp = c(2, 0.7, 0), family = 'Arial', xpd = NA) 
	num.randomisations <- ifelse(is.null(num.randomisations), 
	                             length(randomisation.ids),
	                             num.randomisations)
	# This odd way of getting a sequence is to make sure we get the right outcome if num.randomisations == 0
	r.columns <- randomisation.ids[seq(1, num.randomisations, length.out = num.randomisations)]
	h12.columns <- c('H12', r.columns)
	max.y <- max(c(max(h12.table[, ..h12.columns]), 0.05))
	min.y <- min(h12.table[, ..h12.columns])
	# Create the empty plot.
	plot(c(start.pos, end.pos), c(min.y, max.y), type = 'n', 
	     bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', 
	     main = plot.title, cex.main = 1.2)
	
	# Add randomised data
	sapply(r.columns, function(x) h12.table[, lines(midpoint, get(x), col = colours['randomisations'], lwd = 0.4)])
	# Add true data
	h12.table[, lines(midpoint, H12, col = colours['h12'], lwd = 1.2)]
	# Add y axis
	h12.step.size <- ifelse(max.y > 0.2, 0.1, 0.05)
	axis(2, at = seq(0, max.y, h12.step.size))
	mtext('H12 difference', 2, 2, cex = 0.8)
	# Add peaks that pass filtering. Colour according to whether they have significant p-value
	if (length(p.thresh) == 1)
		p.thresh <- c(p.thresh, -1)
	p.colours <- c('orchid3', 'green', 'blue')
	filter.pass <- h12.table[[filter.name]]
	h12.table[filter.pass, points(midpoint, H12, 
	                              pch = 21,
	                              bg = p.colours[(pval < p.thresh[1]) + (pval < p.thresh[2]) + 1], 
	                              col = colours['randomisations'], 
	                              cex = 1.1, lwd = .5)
	]
	
	# Now add the beadplot
	par(mar = c(2,4,0,1), mgp = c(1, 0.3, 0), tcl = -0.15) 
	plot(c(start.pos, end.pos), c(-2, 0), type = 'n', 
	     bty = 'n', xaxt = 'n', yaxt = 'n', 
	     xlab = paste('Position on', chrom), ylab = '')
    axis(1, at = signif(seq(signif(start.pos, 4), signif(end.pos, 4), length.out = 4), 5))
	draw.gene.model(
		start.pos, 
		end.pos, 
		gene.gff, 
		exon.gff, 
		y = 0, 
		text.cex = 0.5,
		gene.thickness.fraction = 15
	)
	#mtext(paste('Position on', chrom), 1, outer = T)
	
	if (!is.null(filename))
		dev.off()
}

# Create all the plots
lapply(study.pops, 
       function(pop) 
           h12.tables[[pop]]$diff[,
               plot.h12.diff.region(.SD, 
                                    filename = paste(pop, '_region_', .BY$peak.region, '.png', sep = ''),
                                    chrom = .BY$chromosome,
                                    p.thresh = ..p.threshold,
                                    plot.title = paste(pop, 'region', .BY$peak.region, sep = '_')
               ), by = c('peak.region', 'chromosome')
           ]
)

# And we're done.



