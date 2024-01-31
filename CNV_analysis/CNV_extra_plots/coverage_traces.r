library(magrittr)

load('../Ag1000G_CNV_data/v3.7_1246-VO-TZ-KABULA-VMF00185/target_regions_analysis/target_regions_analysis.Rdata')

lighten.col <- function(color, lightness = 1, alpha = alpha){
	col.rgb <- col2rgb(color)/255
	rgb(t(1-(1-col.rgb)*lightness), alpha = alpha)
}

elegant.plot <- function(sample.names, breakpoints, hmm.data, region.coords, gene.coords, 
	                     Dup.name, smoothing = 1, maxy = 12, y.label = '', 
	                     coverage.col = rgb(0.4, 0.4, 0.4), hmm.col = 'brown'){
	par(mar = c(1, 2, 0, 1), mgp = c(1.5,0.5,0))
	plot(region.coords, c(-5,maxy), type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = y.label, cex.axis = 0.8, cex.lab = 0.9, xaxs = 'i')
	region.indices <- which(hmm.data[[1]]$Position > region.coords[1] & hmm.data[[1]]$Position < region.coords[2])
	tickpoints <- c(ceiling(region.coords[1]/10000)*10000, floor(region.coords[2]/10000)*10000)
	#axis(1, at = tickpoints, labels = tickpoints/1000000, tick = F, mgp = c(0,-0.2,3), cex.axis = 0.7, hadj = 0.85)
#	lines(rep(tickpoints[1], 2), c(0, -0.3))
#	lines(rep(tickpoints[2], 2), c(0, -0.3))
	axis(2, at = seq(2, 12, 2))
	abline(h = -0.5, lwd = 2)
	gene.col <- rgb(0.68, 0.84, 0.9, 0.5)
	for (i in seq_len(nrow(gene.coords))){
		rect(gene.coords[i,'start'], -0.5, gene.coords[i, 'end'], maxy, col = gene.col, border = gene.col)
		rect(gene.coords[i,'start'], -0.7, gene.coords[i, 'end'], -0.3, col = 'dodgerblue3', border = 'dodgerblue3')
		gene.name <- rownames(gene.coords)[i]
		if (!grepl('AGAP', gene.name)){
			middle.pos <- mean(c(gene.coords[i, 'start'], gene.coords[i,'end']))
			text(middle.pos, -1, gene.name, srt = 45, adj = 1, font = 2)
		}
	}
	rect(breakpoints[1], -4.7, breakpoints[2], -4.3, col = 'purple', border = 'purple', xpd = NA)
	text(breakpoints[1], -4.5, Dup.name, adj = 1.2, col = 'purple', font = 2)
	for (s in sample.names){
		# If needed, create a table of smoothed data here
		if (smoothing > 1){
			if (smoothing > length(region.indices))
				stop('Fail. Smoothing window size is larger than the number of points to smooth over.')
			these.counts <- hmm.data[[s]][0,]
			for (j in 1:(length(region.indices) - smoothing + 1))
				these.counts[j,] <- apply(hmm.data[[s]][region.indices[j:(j + smoothing - 1)], ], 2, mean)
		}
		else{
			these.counts <- hmm.data[[s]][region.indices, ]
		}
		#lines(these.counts$Position, these.counts$Normalised_coverage, lwd = 2, col = coverage.col)
		points(these.counts$Position, these.counts$Normalised_coverage, col = lighten.col(coverage.col, 0.5), 
		       bg = lighten.col(coverage.col, alpha = 0.3), pch = 21, cex = 1)
		lines(these.counts$Position, these.counts$CNV, lwd = 2, col = hmm.col)
	}
}

Dup0.samples <- rownames(read.based.cnvs$coeaexg)[read.based.cnvs$coeaexg[, 'Dup0']]
Dup1.samples <- rownames(read.based.cnvs$coeaexg)[read.based.cnvs$coeaexg[, 'Dup1']] %>%
                setdiff(high.var.samples) %>%
                intersect(Dup0.samples)
Dup1.breakpoint <- c(known.cnvs$coeaexg$Dup1$BP[1,'pos'], 37295500)
# Don't know why there is a weird aberration in the gene coordinates
gene.coords$coeaexg['AGAP006729', 'Chrom'] <- '2L'
gene.coords$coeaexg['AGAP006729', 'start'] <- 37298660
gene.coords$coeaexg['AGAP006729', 'end'] <- 37299911
coeaexg.region <- c(min(compact.hmm$coeaexg[[1]]$Position), max(compact.hmm$coeaexg[[1]]$Position))
elegant.plot(Dup1.samples[11], Dup1.breakpoint, compact.hmm$coeaexg, coeaexg.region, 
             gene.coords$coeaexg, 'Coeaexg_Dup1')

# Now load the results of the West African data to plot that CNV. For reference, the example 
# sample that I used was VBS19661-5563STDY7800739
compact.hmm.wa <- list(example.wa.sample = read.table('KB_example_sample_coverage.csv', sep = '\t', header = T))

pdf("example_plots.pdf", width = 6.5, height = 4)
par(mfrow = c(2,1), cex = 0.5, oma = c(0,2,0,0))
elegant.plot(Dup1.samples[11], Dup1.breakpoint, compact.hmm$coeaexg, coeaexg.region, 
             gene.coords$coeaexg, 'Coeaexg_Dup1')
elegant.plot('example.wa.sample', c(0,0), compact.hmm.wa, coeaexg.region, gene.coords$coeaexg, '')
mtext('Normalised coverage / copy number state', side = 2, outer = T, cex = 0.7, line = 0.8, font = 2)
dev.off()

png("example_plots.png", width = 6.5, height = 4, units = 'in', res = 300)
par(mfrow = c(2,1), cex = 0.5)
elegant.plot(Dup1.samples[11], Dup1.breakpoint, compact.hmm$coeaexg, coeaexg.region, 
             gene.coords$coeaexg, 'Coeaexg_Dup1')
elegant.plot('example.wa.sample', c(0,0), compact.hmm.wa, coeaexg.region, gene.coords$coeaexg, '')
dev.off()

