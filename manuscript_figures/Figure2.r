# This script will create a summary table of all the genomic regions implicated in resistance based
# on each of the genome-wide anlysis methods. 
library(data.table)
library(stringr)
library(ape)

study.pops <- setNames(nm = c('Moshi_arabiensis_Delta',
                              'Muleba_arabiensis_Delta'))

# We also need a table of region extents. 
gwas.regions <- fread('../misc_scripts/GAARD_SNP/summary_figures/classical_analysis_snp_clump_regions_tanzania.csv') %>%
               .[order(chrom, start)] %>%
			   .[, pos := (start + end)/2] %>%
	           split(by = 'sample.set') %>%
			   .[study.pops]

# Next, the Fst. The results in the haplpotypes folder are derived from the significant Fst regions.
fst.regions <- fread('../haplotypes/haplotype_significance_tests_tanzania.csv') %>%
                     .[, population := factor(population, levels = ..study.pops)] %>%
                     .[!is.na(population)] %>%
                     .[, c('chrom', 'start', 'end') := tstrsplit(window, ':|-')] %>%
                     .[(direction == 'alive') & (logregP < 0.05)] %>%
                     .[, .(population, chrom, start = as.numeric(start), end = as.numeric(end))] %>%
                     unique() %>%
					 split(by = 'population')

# Get the H12 implicated genes. 
h12.p.thresh <- 0.01
h12.regions.fn <- '../GAARD_east/randomisations/H12/h12_filtered_windows_tanzania.RDS'
h12.regions <- readRDS(h12.regions.fn) %>%
               lapply(function(D) D$diff[is.peak == T & pval < h12.p.thresh, 
                                         .(chrom = chromosome, pos = midpoint, start = startpoint, end = endpoint)]
               )

# Get the PBS implicated genes. 
pbs.p.thresh <- 0.01
pbs.regions.fn <- '../GAARD_east/randomisations/PBS/pbs_filtered_windows_tanzania.RDS'
pbs.regions <- readRDS(pbs.regions.fn) %>%
               lapply(function(D) D[is.peak == T & pval < pbs.p.thresh, 
                                    .(chrom = chromosome, pos = midpoint, start = startpoint, end = endpoint)]
               )


# Now we want to plot a figure that summarises the genomic regions implicated by each method. 
chrom.sizes <- c('2R' = 61545105, '2L' = 49364325, '3R' = 53200684, '3L' = 41963435, 'X' = 24393108)

plot.regions.on.chromosome <- function(regions, 
	                                   chrom.sizes, gaps = 5e6,
	                                   show.chrom.names = F, label.cex = 1,
	                                   chrom.col = NULL, chrom.cex = 1.4,
	                                   chrom.thick = 1,
	                                   centro.thick = 0.2,
	                                   chrom.offset = 0, rect.col = 'blue'){
	regions <- copy(regions)
	ce <- cumsum(chrom.sizes + c(0, 0, gaps, 0, gaps))
	cs <- ce - chrom.sizes
	plot(c(cs[1], ce[5]), c(-6.5,1.3), xlim = c(cs[1] + gaps/2, ce[5] - gaps/2), type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
	
	# Add the regions
	if (nrow(regions) > 0){
		regions[, c('genome.start', 'genome.end') := .SD + cs[.BY$chrom], 
				  .SDcols = c('start', 'end'), 
				  by = 'chrom'
		]
		# Give every region a minimum width so that it's visible on the plot
		min.region.size <- 1e6
		regions[, region.size.deficit := ..min.region.size - (end - start)]
		regions[region.size.deficit > 0, ':='(genome.start = genome.start - region.size.deficit/2,
											  genome.end = genome.end + region.size.deficit/2)]
		
		regions[, rect(genome.start, -chrom.thick, genome.end, chrom.thick, 
					   border = NA, col = rect.col, lwd = 5)]
	}
	
	# Draw the outline of the chromosomes
	chrom.col <- if (!is.null(chrom.col)) chrom.col else c('2R' = 'black', '2L' = 'black', '3R' = 'black', '3L' = 'black', 'X' = 'black')
	chrom.y <- -2.5 - chrom.offset
	lines(c(ce['2R'], ce['2R'] - gaps/2, cs['2R'], cs['2R'], ce['2R'] - gaps/2, ce['2R']), 
		  c(-centro.thick, -chrom.thick, -chrom.thick, chrom.thick, chrom.thick, centro.thick), lwd = 2, col = chrom.col['2R'])
	lines(c(cs['2L'], cs['2L'] + gaps/2, ce['2L'], ce['2L'], cs['2L'] + gaps/2, cs['2L']), 
		  c(centro.thick, chrom.thick, chrom.thick, -chrom.thick, -chrom.thick, -centro.thick), lwd = 2, col = chrom.col['2L'])
	lines(c(ce['3R'], ce['3R'] - gaps/2, cs['3R'], cs['3R'], ce['3R'] - gaps/2, ce['3R']), 
		  c(-centro.thick, -chrom.thick, -chrom.thick, chrom.thick, chrom.thick, centro.thick), lwd = 2, col = chrom.col['3R'])
	lines(c(cs['3L'], cs['3L'] + gaps/2, ce['3L'], ce['3L'], cs['3L'] + gaps/2, cs['3L']), 
		  c(centro.thick, chrom.thick, chrom.thick, -chrom.thick, -chrom.thick, -centro.thick), lwd = 2, col = chrom.col['3L'])
	lines(c(cs['X'], cs['X'], ce['X'] - gaps/2, ce['X'], ce['X'], ce['X'] - gaps/2, cs['X']), 
		  c(-chrom.thick, chrom.thick, chrom.thick, centro.thick, -centro.thick, -chrom.thick, -chrom.thick), lwd = 2, col = chrom.col['X'])
	if (show.chrom.names){
		text((cs['2R'] + ce['2R'])/2, chrom.y, '2R', adj = 0.5, xpd = NA, cex = chrom.cex)
		text((cs['2L'] + ce['2L'])/2, chrom.y, '2L', adj = 0.5, xpd = NA, cex = chrom.cex)
		text((cs['3R'] + ce['3R'])/2, chrom.y, '3R', adj = 0.5, xpd = NA, cex = chrom.cex)
		text((cs['3L'] + ce['3L'])/2, chrom.y, '3L', adj = 0.5, xpd = NA, cex = chrom.cex)
		text((cs['X'] + ce['X'])/2, chrom.y, 'X', adj = 0.5, xpd = NA, cex = chrom.cex)
	}
	if (!is.null(regions[['label.offset']]))
		regions[is.na(label.offset), label.offset := 0]
	else
		regions[, label.offset := 0]
	# Add any required labels. We can't use regions$label here, because if the the "label.below" column is
	# present, then "label" will be interpreted as that, and it will think the column is present. Wow. 
	if (!is.null(regions[['label']]))
		regions[, text(genome.end-1e6 + label.offset, 1.5, label, col = rect.col, srt = 35, xpd = NA, adj = 0, cex = label.cex, font = 2)]
	# Add any required labels
	if (!is.null(regions[['label.below']]))
		regions[, text(genome.end + label.offset, -1.5, label.below, col = rect.col, srt = 35, xpd = NA, adj = 1, cex = label.cex, font = 2)]
}

plot.implicated.regions <- function(gwas.regions,
                                    fst.regions,
                                    h12.regions,
                                    pbs.regions,
									study.pops,
									gene.cex = 1,
									label.cex = 1,
									chrom.thick = 0.7,
									centro.thick = 0.2){
	num.pops <- length(study.pops)
	layout(matrix(1:(num.pops*6 - 1)), heights = rep(c(2,6,6,6,6,1.5), num.pops)[-(num.pops*6)])
	par(mar = c(0,3.5,0.6,0), oma = c(0,0,0.1,0), las = 1)
	for (pop in study.pops){
		plot(0,0,type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
		plot.title <- sub('_.*_', ' ', pop) %>%
		              sub('Delta', 'deltamethrin', .)
		mtext(plot.title, line = -0.5, font = 2, cex = 1)
		plot.regions.on.chromosome(gwas.regions[[pop]], chrom.sizes, rect.col = 'purple', label.cex = label.cex, 
		                           chrom.thick = chrom.thick, centro.thick = centro.thick)
		mtext('GWAS', line = -0.5, col = 'purple', font = 2, side = 2, padj = -1)
		plot.regions.on.chromosome(fst.regions[[pop]], chrom.sizes, rect.col = 'orangered3', label.cex = label.cex,
		                           chrom.thick = chrom.thick, centro.thick = centro.thick)
		mtext(expression(F[ST]), line = -0.5, col = 'orangered3', font = 2, side = 2, padj = -1)
		plot.regions.on.chromosome(h12.regions[[gsub('_', '.', pop)]], chrom.sizes, rect.col = 'darkgreen', label.cex = label.cex,
		                           chrom.thick = chrom.thick, centro.thick = centro.thick)
		mtext(expression(paste(Delta, 'H'[12])), line = -0.5, col = 'darkgreen', font = 2, side = 2, padj = -1)
		if (pop != study.pops[num.pops]) {
			plot.regions.on.chromosome(pbs.regions[[gsub('_', '.', pop)]], chrom.sizes, rect.col = 'blue', label.cex = label.cex, 
									   chrom.cex = 1.5, chrom.offset = 2, chrom.thick = chrom.thick, 
									   centro.thick = centro.thick, show.chrom.names = F)
			mtext('PBS', line = -0.5, col = 'blue', font = 2, side = 2, padj = -1)
			plot(0,0,type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
		}
		else {
			plot.regions.on.chromosome(pbs.regions[[gsub('_', '.', pop)]], chrom.sizes, rect.col = 'blue', label.cex = label.cex, 
									   chrom.cex = 1.5, chrom.offset = 2, chrom.thick = chrom.thick, 
									   centro.thick = centro.thick, show.chrom.names = T)
			mtext('PBS', line = -0.5, col = 'blue', font = 2, side = 2, padj = -1)
		}
	}
}


# Manually add some labels to the regions
fst.regions[['Muleba_arabiensis_Delta']][chrom == '2R' & start == 28322651, label := 'Cyp6aa1\n(-120Kb)\n']
# Change the position of that region to fit better on the plot
fst.regions[['Muleba_arabiensis_Delta']][chrom == '2R' & start == 56192837, label := 'Gstd3\n(+700Kb)\n']
#
pbs.regions[['Moshi.arabiensis.Delta']][chrom == '2R' & start == 28560507, label.below := '\nCyp6aa1\n(+80Kb)']
pbs.regions[['Moshi.arabiensis.Delta']][chrom == '2L' & start == 3157513, label.below := '\nVgsc\n(+700Kb)']
pbs.regions[['Muleba.arabiensis.Delta']][chrom == '2R' & start == 28327171, label := 'Cyp6aa1\n(-110Kb)\n']
pbs.regions[['Muleba.arabiensis.Delta']][chrom == '2R' & start == 56707016, label := 'Gstd3\n(+1.2Mb)\n']
#
h12.regions[['Moshi.arabiensis.Delta']][chrom == '2R' & start == 28536270, label := 'Cyp6aa1']
h12.regions[['Moshi.arabiensis.Delta']][chrom == '2R' & start == 40837427, label.below := '\nKeap1\n(-65Kb)']
h12.regions[['Moshi.arabiensis.Delta']][chrom == '2L' & start == 28491089, label := 'Coeae1f\n(-6Kb)\n']
h12.regions[['Moshi.arabiensis.Delta']][chrom == '2L' & start == 28491089, label.offset := -2e6]
h12.regions[['Moshi.arabiensis.Delta']][chrom == '2L' & start == 36821514, label.below := '\nCoeae2g-7g\n(-440Kb)']
h12.regions[['Muleba.arabiensis.Delta']][chrom == '2R' & start == 28416418, label := 'Cyp6aa1']
h12.regions[['Muleba.arabiensis.Delta']][chrom == '2R' & start == 20316806, label.below := '\nCyp4k2\n(-700Kb)']
h12.regions[['Muleba.arabiensis.Delta']][chrom == '2R' & start == 20316806, label.offset := -1e6]
#
#gwas.regions[['Muleba_arabiensis_Delta']][chrom == '2R' & start == 28650000, label.below := 'Cyp6aa1 (+170Kb)']
gwas.regions[['Muleba_arabiensis_Delta']][chrom == '2R' & start == 28650000, label := 'Cyp6aa1\n(+170Kb)\n']
gwas.regions[['Muleba_arabiensis_Delta']][chrom == '3R' & start == 9150000, label.below := '\nCyp4h16\n(-600Kb)']
#gwas.regions[['Moshi_arabiensis_Delta']][chrom == '3R' & start == 5050000, label.below := 'Cyp12f1 (+800Kb)']
gwas.regions[['Moshi_arabiensis_Delta']][chrom == '3R' & start == 5050000, label.below := '\nCyp12f1\n(+800Kb)']

filename.png <- paste('Figure2.png', sep = '_')
png(filename.png, width = 7, height = 5.6, units = 'in', res = 300)
plot.implicated.regions(gwas.regions, 
						fst.regions,
						h12.regions,
						pbs.regions,
						study.pops,
						gene.cex = 1,
						chrom.thick = 0.6)
dev.off()

filename.tiff <- paste('Fig2.tiff', sep = '_')
tiff(filename.tiff, width = 7, height = 5.6, units = 'in', res = 300)
plot.implicated.regions(gwas.regions, 
						fst.regions,
						h12.regions,
						pbs.regions,
						study.pops,
						gene.cex = 1,
						chrom.thick = 0.6)
dev.off()

filename.svg <- paste('Figure2.svg', sep = '_')
svg(filename.svg, width = 7, height = 5.6)
plot.implicated.regions(gwas.regions, 
						fst.regions,
						h12.regions,
						pbs.regions,
						study.pops,
						gene.cex = 1,
						chrom.thick = 0.6)
dev.off()



