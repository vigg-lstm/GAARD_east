library(data.table)

snp.frequencies.tanz <- fread('target_site_genotype_frequencies_tanzania.csv', key = 'Locus')

snp.freq.flipped.tanz <- t(snp.frequencies.tanz[, -1])
snp.names.tanz <- snp.frequencies.tanz[[1]]

# All SNPs can have the mutant nucleotide part of the name removed
snp.names.tanz <- sub('_.*', '', snp.names.tanz)
colnames(snp.freq.flipped.tanz) <- snp.names.tanz

source('../shared_functions/R_plotting.r')

# A function that will round to x significant digits, unless the number is small, 
# then just take one
signif2 <- function(x, digits, threshold = 0.1){
	small.x <- x < threshold 
	x[small.x] <- signif(x[small.x], 1)
	x[!small.x] <- signif(x[!small.x], digits)
	x
}

plot.table <- function(freq.table, species.colours = c(coluzzii = 'coral3', gambiae = 'dodgerblue3', arabiensis = 'limegreen'), 
                       text.cell.cex = 0.65, pop.cex = 0.5, snp.cex = 0.5, ...){
	if ((length(names(species.colours)) == 0) | !any(names(species.colours) %in% c('arabiensis', 'gambiae', 'coluzzii', 'intermediate')))
		names(species.colours) <- c('arabiensis', 'gambiae', 'coluzzii', 'intermediate')
	par(mar = c(0,3,3,0), ...)
	plot(c(0,ncol(freq.table)), -c(0,nrow(freq.table)), type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
	x.matrix <- matrix(1:ncol(freq.table), nrow(freq.table), ncol(freq.table), byrow = T)
	y.matrix <- matrix(1:nrow(freq.table), nrow(freq.table), ncol(freq.table))
	text.col.matrix <- c('black', 'white')[(freq.table > 0.5) + 1]
	text.cell.matrix <- signif2(freq.table, 2)
	col.matrix <- matrix(rgb(0.8, 0.8, 0.8), nrow(freq.table), ncol(freq.table))
	for (i in 1:nrow(freq.table)){
		this.row <- freq.table[i, ]
		if (any(!is.na(this.row))){
			this.species <- sub('^.*\\.', '', rownames(freq.table)[i])
			col.matrix[i, !is.na(this.row)] <- lighten.col(species.colours[this.species], sqrt(this.row[!is.na(this.row)]))
	
		}
	}
	rect(x.matrix - 1, 1 - y.matrix, x.matrix, - y.matrix, col = col.matrix, border = 'white', lwd = text.cell.cex)
	text(x.matrix - 0.5, 0.5 - y.matrix, text.cell.matrix, cex = text.cell.cex, col = text.col.matrix, font = 2)
	text(x.matrix[1,] - 0.5, (par('usr')[4] - par('usr')[3])/100, srt = 45, adj = c(0,0), sub('\\.', ' ', colnames(freq.table)), xpd = NA, cex = snp.cex)
	country.names <- sub('\\.[^.]*$', '', rownames(freq.table))
	species <- sub('.*\\.', '', rownames(freq.table))
	rect(x.matrix[1,1]-1, 1 - y.matrix, par('usr')[1] - (par('usr')[2] - par('usr')[1])/5, -y.matrix, col = sapply(species.colours[species], lighten.col, 0.4), border = NA, xpd = NA)
	text(x.matrix[1,1]-1, 0.5 - y.matrix[,1], paste(country.names, '  \n', species, '  ', sep = ''), adj = 1, cex = pop.cex, xpd = NA)
}


png('target_snp_frequencies.png', width = 2, height = 0.8, units = 'in', res = 300)
par(family = 'Arial')
plot.table(snp.freq.flipped.tanz, 
           text.cell.cex = 0.35,
           pop.cex = 0.3,
           snp.cex = 0.3,
           mai = c(0,0.2,0.22,0)
)
dev.off()

