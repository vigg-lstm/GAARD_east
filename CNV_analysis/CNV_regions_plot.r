library(stringr)
library(data.table)
library(ape)

# Load the gff
gff.path <- '../data/VectorBase-57_AgambiaePEST.gff'
gff <- as.data.table(read.gff(gff.path, GFF3 = T))
gff[, seqid := sub('AgamP4_', '', seqid)]

# Get the list of gene present in the Cyp6 region
cyp6.plotting.region <- c(28465000, 28572600)
cyp6.gff <- gff[seqid == '2R' & 
				start < cyp6.plotting.region[2] & 
				end > cyp6.plotting.region[1] &
				type == 'protein_coding_gene', ] %>%
			.[order(start)]
cyp6.gff[, gene.id := unique(str_extract(attributes, 'AGAP\\d{6}'))]
cyp6.gff[, gene.name := str_extract(attributes, '(?<=Name=)[^;]+') %>%
						str_to_title]

# Only keep the genes that have names (this kicks out anything outside of the Cyp6aa / Cyp6p cluster
cluster.cyp6.gff <- cyp6.gff[!is.na(gene.name)]

# The new CNV ranges:
new.cyp6.cnvs = list(Dup31 = c(28480476, 28492867),
                     Dup32 = c(28485349, 28494719),
                     Dup33 = c(28479842, 28549886),
                     Dup34 = c(28480284, 28483348),
                     Dup35 = c(28472668, 28504018),
                     Dup36 = c(28475403, 28572600),
                     Dup37 = c(28477192, 28489472))

gene.line.pos <- -0.3
gene.text.offset <- 0.3
gene.cex = 0.7
dup.col <- rgb(0.79, 0.7, 0.84)
dup.border.col <- 'orchid'
#dup.col <- 'pink'
#dup.border.col <- 'red'
dup.text.col <- 'black'
dup.height <- 0.6
dup.cex <- 0.75
dup.text.adjustment <- 0.4

# Right, now we want to plot the cyp6 region, along with the various duplication types
#png('new_cyp6_cnvs.png', width = 12.5, height = 3.5, units = 'in', res = 400)
pdf('new_cyp6_cnvs.pdf', width = 12.5, height = 3.5)
par(mar = c(4,1.5,1,1), mgp = c(2, 0.7, 0))
plot(cyp6.plotting.region, c(2.2,-length(new.cyp6.cnvs)-0.5), type = 'n', xlab = 'Position on Chromosome 2R', ylab = '', xaxt = 'n', yaxt = 'n', cex.lab = 1.1) 
lines(cyp6.plotting.region, rep(gene.line.pos, 2), lwd = 2)
# Show the position of the genes of interest
for (i in seq_len(nrow(cluster.cyp6.gff))){
	gene.name <- cluster.cyp6.gff[i, gene.name]
	if (gene.name == 'Cyp6aa1'){
		gene.col <- 'orange'
		bar.col <- rgb(1,1,0.5)
	}
	else if (gene.name == 'Cyp6p3'){
		gene.col <- 'seagreen'
		bar.col <- 'lightgreen'
	}	
	else {
		gene.col <- 'black'
		if ((i %% 2) == 0)
			bar.col <- 'grey80'
		else
			bar.col <- 'grey60'
	}
	rect(cluster.cyp6.gff[i, start], -100, cluster.cyp6.gff[i, end], 100, col = bar.col, border = bar.col)
	rect(cluster.cyp6.gff[i, start], gene.line.pos-0.15, cluster.cyp6.gff[i, end], gene.line.pos+0.15, col = gene.col, border = gene.col)
	text(cluster.cyp6.gff[i, mean(c(start, end))], gene.line.pos+gene.text.offset, gene.name, srt=90, adj = 0, col = 'black', cex = gene.cex, font = 2)
}
# Draw the axis (also redraw the top axis to cover encroaching plot).
axis(1, at = seq(round(cyp6.plotting.region[1], -4), round(cyp6.plotting.region[2], -4), 2e4), cex.axis = 0.8) 
axis(3, labels = F, cex.axis = 0.8, tcl = 0) 
# plot the duplications of interest
for (i in seq_along(new.cyp6.cnvs)){
	cnv.name <- names(new.cyp6.cnvs[i])
	cnv.range <- new.cyp6.cnvs[[i]] 
	rect(cnv.range[1], -(i + dup.height), cnv.range[2], -i, col = dup.col, border = dup.border.col, lwd = 0.3)
	text(cnv.range[2], -(i + dup.text.adjustment), cnv.name, cex = dup.cex, adj = c(1,0.1), font = 2, col = dup.text.col)
}
dev.off()

# Now Coeaexf
coeaexf.plotting.region <- c(28530000, 28572000)
coeaexf.gff <- gff[seqid == '2L' & 
				start < coeaexf.plotting.region[2] & 
				end > coeaexf.plotting.region[1] &
				type == 'protein_coding_gene', ] %>%
			.[order(start)]
coeaexf.gff[, gene.id := unique(str_extract(attributes, 'AGAP\\d{6}'))]
coeaexf.gff[, gene.name := str_extract(attributes, '(?<=Name=)[^;]+') %>%
						str_to_title]
# Coeae1f is not named in the gff
coeaexf.gff[gene.id == 'AGAP006227', 'gene.name'] <- 'Coeae1f'

new.coeaexf.cnvs = list(Dup2 = c(28535653, 28571586),
                        Dup3 = c(28533512, 28550548),
                        Dup4 = c(28533499, 28557804),
                        Dup5 = c(28541709, 28554102),
                        Dup6 = c(28548424, 28552841),
                        Dup7 = c(28546647, 28552722))

dup.col <- 'pink'
dup.border.col <- 'red'

pdf('new_coef_cnvs.pdf', width = 12.5, height = 3.5)
par(mar = c(4,1.5,1,1), mgp = c(2, 0.7, 0))
plot(coeaexf.plotting.region, c(2.5,-length(new.coeaexf.cnvs)-0.5), type = 'n', xlab = 'Position on Chromosome 2L', ylab = '', xaxt = 'n', yaxt = 'n', cex.lab = 1.1) 
lines(coeaexf.plotting.region, rep(gene.line.pos, 2), lwd = 2)
# Show the position of the genes of interest
for (i in seq_len(nrow(coeaexf.gff))){
	gene.name <- sub('_', '\n', coeaexf.gff[i, gene.name])
	if (is.na(gene.name))
		gene.name <- coeaexf.gff[i, gene.id]
	if (gene.name %in% c('Coeae1f', 'Coeae2f')){
		gene.col <- 'blue'
		bar.col <- rgb(0.8,0.8,1)
	}
	else {
		gene.col <- 'black'
		if ((i %% 2) == 0)
			bar.col <- 'grey80'
		else
			bar.col <- 'grey60'
	}
	rect(coeaexf.gff[i, start], -100, coeaexf.gff[i, end], 100, col = bar.col, border = bar.col)
	rect(coeaexf.gff[i, start], gene.line.pos-0.15, coeaexf.gff[i, end], gene.line.pos+0.15, col = gene.col, border = gene.col)
	text(coeaexf.gff[i, mean(c(start, end))], gene.line.pos+gene.text.offset, gene.name, srt=90, adj = 0, col = 'black', cex = gene.cex, font = 2)
}
# Draw the axis (also redraw the top axis to cover encroaching plot).
axis(1, at = seq(round(coeaexf.plotting.region[1], -4), round(coeaexf.plotting.region[2], -4), 1e4), cex.axis = 0.8) 
axis(3, labels = F, cex.axis = 0.8, tcl = 0) 
# plot the duplications of interest
for (i in seq_along(new.coeaexf.cnvs)){
	cnv.name <- names(new.coeaexf.cnvs[i])
	cnv.range <- new.coeaexf.cnvs[[i]] 
	rect(cnv.range[1], -(i + dup.height), cnv.range[2], -i, col = dup.col, border = dup.border.col, lwd = 0.3)
	text(cnv.range[2], -(i + dup.text.adjustment), cnv.name, cex = dup.cex, adj = c(1,0.1), font = 2, col = dup.text.col)
}
dev.off()


# Now Coeaexg
coeaexg.plotting.region <- c(37280000, 37300000)
coeaexg.gff <- gff[seqid == '2L' & 
				start < coeaexg.plotting.region[2] & 
				end > coeaexg.plotting.region[1] &
				type == 'protein_coding_gene', ] %>%
			.[order(start)]
coeaexg.gff[, gene.id := unique(str_extract(attributes, 'AGAP\\d{6}'))]
coeaexg.gff[, gene.name := str_extract(attributes, '(?<=Name=)[^;]+') %>%
						str_to_title]

# Coeae7g is annotated in panoptes, but doesn't seem to be here. Add the annotation
coeaexg.gff[gene.id == 'AGAP006728', gene.name := 'Coeae7g']
new.coeaexg.cnvs = list('Dup1' = c(37282078, 37295500),
                        'Dup2' = c(37294010, 37299365))


dup.col <- 'salmon'
dup.border.col <- 'brown'

pdf('new_coeg_cnvs.pdf', width = 12.5, height = 2.5)
par(mar = c(4,1.5,1,1), mgp = c(2, 0.7, 0)) 
plot(coeaexg.plotting.region, c(2.5,-length(new.coeaexg.cnvs)-0.5), type = 'n', xlab = 'Position on Chromosome 2L', ylab = '', xaxt = 'n', yaxt = 'n', cex.lab = 1.1) 
lines(coeaexg.plotting.region, rep(gene.line.pos, 2), lwd = 2)
# Show the position of the genes of interest
for (i in seq_len(nrow(coeaexg.gff))){
	gene.name <- sub('_', '\n', coeaexg.gff[i, gene.name])
	if (is.na(gene.name))
		gene.name <- coeaexg.gff[i, gene.id]
	gene.col <- 'black'
	if ((i %% 2) == 0)
		bar.col <- 'grey80'
	else
		bar.col <- 'grey60'
	rect(coeaexg.gff[i, start], -100, coeaexg.gff[i, end], 100, col = bar.col, border = bar.col)
	rect(coeaexg.gff[i, start], gene.line.pos-0.15, coeaexg.gff[i, end], gene.line.pos+0.15, col = gene.col, border = gene.col)
	text(coeaexg.gff[i, mean(c(start, end))], gene.line.pos+gene.text.offset, gene.name, srt=90, adj = 0, col = 'black', cex = gene.cex, font = 2)
}
# Draw the axis (also redraw the top axis to cover encroaching plot).
axis(1, at = seq(round(coeaexg.plotting.region[1], -4), round(coeaexg.plotting.region[2], -4), 1e4), cex.axis = 0.8) 
axis(3, labels = F, cex.axis = 0.8, tcl = 0) 
# plot the duplications of interest
for (i in seq_along(new.coeaexg.cnvs)){
	cnv.name <- names(new.coeaexg.cnvs[i])
	cnv.range <- new.coeaexg.cnvs[[i]] 
	rect(cnv.range[1], -(i + dup.height), cnv.range[2], -i, col = dup.col, border = dup.border.col, lwd = 0.3)
	text(cnv.range[2], -(i + dup.text.adjustment), cnv.name, cex = dup.cex, adj = c(1,0.1), font = 2, col = dup.text.col)
}
dev.off()


