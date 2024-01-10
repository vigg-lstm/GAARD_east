library(stringi)
library(stringr)
library(data.table)
library(magrittr)
library(ape)

# Load all of the H12 tables
h12.filenames <- list.files('.', 'H12_.*.csv', full.names = T)

#study.pops <- setNames(nm = c('Moshi.arabiensis.Delta',
#                              'Moshi.arabiensis.PM',
#                              'Muleba.arabiensis.Delta'))
study.pops <- setNames(nm = c('Moshi',
                              'Muleba'))

# A function to load and combine all data for a given randomisation id
load.and.combine.h12 <- function(pop){
	
	filenames <- grep(pop, h12.filenames, value = T) %>%
	             setNames(., stri_extract_first_regex(., '[23LRX]+(?=\\.csv)'))
	
	combined.data <- names(filenames) %>%
	                 lapply(function(chrom){ 
	                     fread(filenames[chrom], sep = '\t') %>%
	                     .[, chromosome := ..chrom] %>%
					     setnames('h12', 'H12') %>%
					     setcolorder(c('chromosome', 'midpoint'))
				     }) %>%
				     rbindlist()
	
	combined.data
}

cat('Loading H12 data\n')
h12.tables <- lapply(study.pops, function(pop) {cat(pop, '\n'); load.and.combine.h12(pop)})

source('../shared_functions/R_plotting.r')

chrom.sizes <- c('2R' = 61545105, '2L' = 49364325, '3R' = 53200684, '3L' = 41963435, 'X' = 24393108)

# Function to plot the H12 data
plot.h12 <- function(h12.table, 
                     filename = NULL, 
                     plot.title = '', 
                     gaps = 5000000, 
                     h12.col = 'darkgreen'){
	# Create the plot
	if (!missing(filename)){
		file.width = 6.5
		file.height = 2.7
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
	par(mar = c(0,4,1,2), mgp = c(2, 0.7, 0), family = 'Arial', xpd = NA) 
	max.y <- max(c(max(h12.table[, 'H12']), 0.05))
	min.y <- min(h12.table[, 'H12'])
	# Create the empty plot.
	plot(c(cs[1], ce[5]), c(min.y, max.y), xlim = c(cs[1] + gaps/2, ce[5] - gaps/2), 
	     type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', 
	     main = plot.title, cex.main = 1.2)
	
	# Get plotting positions for windows
	h12.table$genome.pos <- h12.table$midpoint + cs[as.character(h12.table$chromosome)]
	# Add the data
	h12.table[, lines(genome.pos, H12, col = h12.col, lwd = 1.2), by = chromosome]
	# Add y axis
	h12.step.size <- ifelse(max.y > 0.2, 0.1, 0.05)
	axis(2, at = seq(0, max.y, h12.step.size))
	mtext('H12 difference', 2, 2, cex = 0.8)
	
	# Now plot all chromosomes with, the position of each of the four detox gene regions and Ace1
	par(mar = c(1,4,0,2), mgp = c(2, 0.7, 0)) 
	add.chromosomes(chrom.sizes, gaps = gaps, gene.cex = 0.6, point.cex = 1, chrom.offset = -1.2, chrom.cex = 1.2)
	
	if (!missing(filename))
		dev.off()
}

for (pop in names(h12.tables))
	plot.h12(h12.tables[[pop]], 
	         filename = paste(pop, 'H12.png', sep = '_'),
	         plot.title = pop)


# Now all the same for H1x
# Load all of the H12 tables
h1x.filenames <- list.files('.', 'H1x_.*.csv', full.names = T)

#pop.pairs <- setNames(nm = c('Moshi.arabiensis.Delta_vs_Moshi.arabiensis.PM',
#                             'Moshi.arabiensis.Delta_vs_Muleba.arabiensis.Delta',
#                             'Moshi.arabiensis.PM_vs_Muleba.arabiensis.Delta'))
pop.pairs <- setNames(nm = c('Moshi_vs_Muleba'))

# A function to load and combine all data for a given randomisation id
load.and.combine.h1x <- function(pop1, pop2){
	
	filenames <- grep(pop1, h1x.filenames, value = T) %>%
	             grep(pop2, ., value = T) %>%
	             setNames(., stri_extract_first_regex(., '[23LRX]+(?=\\.csv)'))
	
	combined.data <- names(filenames) %>%
	                 lapply(function(chrom){ 
	                     fread(filenames[chrom], sep = '\t') %>%
	                     .[, chromosome := ..chrom] %>%
					     setnames('h1x', 'H1x') %>%
					     setcolorder(c('chromosome', 'midpoint'))
				     }) %>%
				     rbindlist()
	
	combined.data
}

cat('Loading data\n')
h1x.tables <- lapply(strsplit(pop.pairs, '_vs_'),
                     function(pops) {print(pops); load.and.combine.h1x(pops[1], pops[2])}
)

# Function to plot the H12 data
plot.h1x <- function(h1x.table, 
                     filename = NULL, 
                     plot.title = '', 
                     gaps = 5000000, 
                     h1x.column = 'H1x',
                     h1x.colour = 'darkgreen',
                     y.label = 'H1X'){
	# Create the plot
	if (!missing(filename)){
		file.width = 6.5
		file.height = 2.7
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
	par(mar = c(0,4,1,2), mgp = c(2, 0.7, 0), family = 'Arial', xpd = NA) 
	max.y <- max(c(max(h1x.table[, ..h1x.column]), 0.05))
	min.y <- min(h1x.table[, ..h1x.column])
	# Create the empty plot.
	plot(c(cs[1], ce[5]), c(min.y, max.y), xlim = c(cs[1] + gaps/2, ce[5] - gaps/2), 
	     type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', 
	     main = plot.title, cex.main = 1.2)
	
	# Get plotting positions for windows
	h1x.table$genome.pos <- h1x.table$midpoint + cs[as.character(h1x.table$chromosome)]
	# Add the data
	h1x.table[, lines(genome.pos, get(h1x.column), col = h1x.colour, lwd = 1.2), by = chromosome]
	# Add y axis
	h1x.step.size <- ifelse(max.y > 0.2, 0.1, 0.05)
	axis(2, at = seq(0, max.y, h1x.step.size))
	mtext(y.label, 2, 2, cex = 0.8)
	
	# Now plot all chromosomes with, the position of each of the four detox gene regions and Ace1
	par(mar = c(1,4,0,2), mgp = c(2, 0.7, 0)) 
	add.chromosomes(chrom.sizes, gaps = gaps, gene.cex = 0.6, point.cex = 1, chrom.offset = -1.2, chrom.cex = 1.2)
	
	if (!missing(filename))
		dev.off()
}

for (pop.pair in names(h1x.tables))
	plot.h1x(h1x.tables[[pop.pair]], 
	         filename = paste(pop.pair, 'H1x.png', sep = '_'),
	         plot.title = pop.pair,
	         h1x.column = 'H1x',
	         h1x.colour = 'purple',
	         y.label = 'H1x')

# Now plot specific regions in H12. Start by manually defining regions of interest. These are
# just indexes that we determined subjectively
regions.of.interest.moshi <- list(region1 = 1075:1150, 
                                  region2 = 1525:1600, 
                                  region3 = 3475:3600, 
                                  region4 = 3690:3850, 
                                  region5 = 4310:4360, 
                                  region6 = 7175:7250)
regions.of.interest.muleba <- list(region1 = 3260:3350, 
                                   region2 = 3700:3825)

# Load the gff
gff.path <- '../data/VectorBase-57_AgambiaePEST.gff'
gff <- as.data.table(read.gff(gff.path, GFF3 = T))
gff[, seqid := sub('AgamP4_', '', seqid)]


plot.h12.region <- function(h12.table, 
                            filename = NULL, 
                            genes.filename = NULL, 
                            chrom = NULL, 
                            plot.title = '', 
	                        h12.col = 'darkgreen'){
	if (is.null(chrom)){
		chrom <- unique(h12.table$chromosome)
		if (length(chrom) != 1)
			stop('Region should include exactly one chromosome')
	}
		
	start.pos <- min(h12.table$midpoint)
	end.pos <- max(h12.table$midpoint)
	
	# Get the list of genes present in the window
	gene.gff <- gff[seqid == chrom & 
	                start < end.pos & 
	                end > start.pos &
	                type == 'protein_coding_gene', ] %>%
	            .[order(start)]
	gene.gff[, gene.id := unique(str_extract(attributes, 'AGAP\\d{6}'))]
	gene.gff[, gene.name := str_extract(attributes, '(?<=Name=)[^;]+') %>%
	                        str_to_title]
	gene.gff$gene.name[is.na(gene.gff$gene.name)] <- gene.gff$gene.id[is.na(gene.gff$gene.name)]
	if (!is.null(genes.filename))
		fwrite(gene.gff[, .(seqid, start, end, gene.id, gene.name, attributes)], genes.filename, sep = '\t')
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
	par(mar = c(0,4,1,1), mgp = c(2, 0.7, 0), family = 'Arial', xpd = NA) 
	# This odd way of getting a sequence is to make sure we get the right outcome if num.randomisations == 0
	max.y <- max(c(max(h12.table[, 'H12']), 0.05))
	min.y <- min(h12.table[, 'H12'])
	# Create the empty plot.
	plot(c(start.pos, end.pos), c(min.y, max.y), type = 'n', 
	     bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', 
	     main = plot.title, cex.main = 1.2)
	
	# Add data
	h12.table[, lines(midpoint, H12, col = h12.col, lwd = 1.2)]
	# Add y axis
	h12.step.size <- ifelse(max.y > 0.2, 0.1, 0.05)
	axis(2, at = seq(0, max.y, h12.step.size))
	mtext('H12', 2, 2, cex = 0.8)
	
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

for (r in regions.of.interest.moshi){
	chrom = unique(h12.tables$Moshi[r, chromosome])
	start.point = h12.tables$Moshi[r[1], midpoint]
	end.point = h12.tables$Moshi[r[length(r)], midpoint]
	fn.root <- paste('Moshi', chrom, start.point, end.point, sep = '_')
	plot.h12.region(h12.tables$Moshi[r], 
	                paste(fn.root, '.png', sep = ''),
	                paste(fn.root, '_genes.csv', sep = '')
	)
}

for (r in regions.of.interest.muleba){
	chrom = unique(h12.tables$Muleba[r, chromosome])
	start.point = floor(h12.tables$Muleba[r[1], midpoint])
	end.point = floor(h12.tables$Muleba[r[length(r)], midpoint])
	fn.root <- paste('Muleba', chrom, start.point, end.point, sep = '_')
	plot.h12.region(h12.tables$Muleba[r], 
	                paste(fn.root, '.png', sep = ''),
	                paste(fn.root, '_genes.csv', sep = '')
	)
}
