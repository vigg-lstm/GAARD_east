library(data.table)
library(magrittr)
library(stringr)
library(glmmTMB)
library(abind)

source('../shared_functions/R_plotting.r')
study.ids.table <- fread('../data/study_ids.csv')
study.ids.table[, study_id := paste('v3.7_', study_id, sep = '')]
meta <- fread('../data/tanzania.samples.meta.csv', key = 'sample_id')
meta.sample.names <- meta$partner_sample_id
phen <- fread('../data/sample_phenotypes_EA.csv', key = 'specimen')[sort(meta.sample.names), ]

# Identify samples for removal (sibs, males and samples to exlude due to unusual
# relatedness)
sib.groups <- fread('../NGSrelate/full_relatedness_tanzania/sib_group_table.csv')
sibs.to.remove <- sib.groups[keep == F, sample.name]
males <- meta[sex_call == 'M', partner_sample_id]
samples.to.exclude = readLines('../NGSrelate/full_relatedness_tanzania/samples_to_exclude.csv')
samples.to.remove <- c(sibs.to.remove, males, samples.to.exclude)

load.target.cnv.table <- function(study.id){
	target.cnv.table <- fread(paste('Ag1000G_CNV_data', study.id, 'target_regions_analysis/focal_region_CNV_table.csv', sep = '/'))
	these.sample.names <- meta[target.cnv.table$V1, partner_sample_id]
	target.cnv.table$sample.id <- these.sample.names
	target.cnv.table$V1 <- NULL
	target.cnv.table
}

target.CNV.table <- lapply(unique(study.ids.table$study_id), load.target.cnv.table) %>%
                    rbindlist() %>%
                    .[!(sample.id %in% samples.to.remove)] %>%
                    setkey(sample.id)

Dup.names <- grep('_D[upel]{2}', colnames(target.CNV.table), value = T)
Dup.families <- sub('_.*', '', Dup.names)
Dup.clusters <- split(Dup.names, Dup.families)

known.Dups <- Dup.names[!grepl('Dup0', Dup.names)]
sample.names <- target.CNV.table$sample.id
good.var.sample.names <- target.CNV.table[High.var.sample == F]$sample.id
# The "Any" category is where any of the known CNVs are present
Dups.any <- c('Ace1_Any', 'Cyp6aap_Any', 'Cyp6mz_Any', 'Gstue_Any', 'Cyp9k1_Any', 'Coeaexf_Any', 'Coeaexg_Any')
for (D in Dups.any){
	this.region <- sub('_Any', '', D)
	these.cnvs <- known.Dups[grepl(this.region, known.Dups)]
	target.CNV.table[[D]] <- apply(target.CNV.table[, ..these.cnvs, drop = F], 1, any)
}

target.CNV.table$phenotype <- as.factor(phen[sample.names, phenotype])

# Get the CNV count per population
population.target.numCNVs <- aggregate(target.CNV.table[, c(Dup.names, Dups.any), with = F], 
                                       phen[sample.names, .(location, species)],
                                       sum
)
# Get the CNV frequency per population
population.target.CNVs <- aggregate(target.CNV.table[, c(Dup.names, Dups.any), with = F], 
                                    phen[sample.names, .(location, species)], 
                                    function(x) sum(x) / length(x)
)
# Get the size of each population-year combination
population.popsize <- interaction(phen[sample.names, .(location, species)]) %>%
                      droplevels %>%
                      table
# Now again using only samples with good variance
population.target.goodvar.numCNVs <- aggregate(target.CNV.table[good.var.sample.names, c(Dup.names, Dups.any), with = F], 
                                               phen[good.var.sample.names, .(location, species)],
                                               sum
)
population.target.goodvar.CNVs <- aggregate(target.CNV.table[good.var.sample.names, c(Dup.names, Dups.any), with = F], 
                                            phen[good.var.sample.names, .(location, species)], 
                                            function(x) sum(x) / length(x)
)
population.goodvar.popsize <- interaction(phen[good.var.sample.names, .(location, species)]) %>%
                      droplevels %>%
                      table

# Function to spin a table of by-population values and fill in missing genes with zeros.
spin.and.fill <- function(x, keyword) {
	out <- as.data.frame(t(x[, c(Dup.names, Dups.any)]))
	out[is.na(out)] <- 0
	colnames(out) <- paste(x$location, x$species, keyword, sep = '.')
	out
}

focal.region.cnv.frequencies <- cbind(spin.and.fill(population.target.numCNVs, 'all'),
                                      spin.and.fill(population.target.CNVs, 'all.prop'),
                                      spin.and.fill(population.target.goodvar.numCNVs, 'goodvar'),
                                      spin.and.fill(population.target.goodvar.CNVs, 'goodvar.prop')
)


# Now get the modal CNVs by gene
load.modal.cnv.table <- function(study.id){
	modal.cnv.table <- fread(paste('Ag1000G_CNV_data', study.id, 'modal_CNVs/modal_copy_number_arabiensis.csv', sep = '/'))
	these.sample.names <- meta[modal.cnv.table$V1, partner_sample_id]
	modal.cnv.table$V1 <- these.sample.names
	colnames(modal.cnv.table)[1] <- 'sample.id'
	modal.cnv.table
}

modal.copy.number.all <- lapply(unique(study.ids.table$study_id), load.modal.cnv.table) %>%
                         rbindlist() %>%
                         .[!(sample.id %in% samples.to.remove)] %>%
                         .[high_variance == F, -c('sex_call', 'high_variance')]

# Since the HMM only goes up to 12 (10 extra copies), we need to recalculate the copy number state using raw
# coverge values for Coeaexg
load.coeaexg.copy.numbers <- function(study.id){
	cn.table <- readRDS(paste('Ag1000G_CNV_data', study.id, 'target_regions_analysis/Coeaexg_cnv_coverage.RDS', sep = '/'))
	these.sample.names <- meta[dimnames(cn.table)[[2]], partner_sample_id]
	dimnames(cn.table)[[2]] <- these.sample.names
	cn.table
}

coeaexg.copy.numbers <- lapply(unique(study.ids.table$study_id), load.coeaexg.copy.numbers) %>%
                        abind(along = 2) %>%
                        .[,!(dimnames(.)[[2]] %in% samples.to.remove),]
coeaexg.median.cn <- round(apply(coeaexg.copy.numbers[, , 'Normalised_coverage'], 2, median)) - 2
modal.copy.number.all[, Coeaexg := ..coeaexg.median.cn[sample.id]]

modal.CNV.table <- copy(modal.copy.number.all) %>%
                   .[, colnames(.)[-1] := lapply(.SD, `>`, 0), .SDcols = colnames(.)[-1]]
# For a smaller table, we only keep genes where a CNV was present. 
no.cnv <- names(which(apply(modal.CNV.table[, -1], 2, sum, na.rm = T) < 1))
modal.CNV.table <- modal.CNV.table[, -c(..no.cnv)]

# Calculate the frequency of CNVs in each population by year
population.modal.CNVs <- aggregate(modal.CNV.table[, -c('sample.id')],
                                   phen[modal.CNV.table$sample.id, .(location, species)], 
                                   function(x) sum(x) / length(x)
)

detox.genes <- c('Ace1',
                 paste('Coeae', 1:2, 'f', sep = ''),
                 'Coeaexg',
                 paste('Cyp6aa', 1:2, sep = ''),
                 paste('Cyp6p', 1:5, sep = ''),
                 paste('Cyp6m', 2:4, sep = ''),
                 paste('Cyp6z', 1:3, sep = ''),
                 paste('Gste', 1:8, sep = ''),
                 'Cyp9k1')

# Get the names of the detox genes
gene.table <- fread('Ag1000G_CNV_data/gene_annotation_fullgenetable.csv', key = 'Gene_stable_ID', check.names = T)
# Coeae1f is absent from the gff, so add it manually
gene.table['AGAP006227', Gene.name := 'COEAE1F']
# We also add a special entry for the Coeaexg cluster, for subsequent functions to work
gene.table <- rbind(gene.table, data.table(Gene_stable_ID = 'Coeaexg', Gene.name = 'COEAEXG', Descriptions = '', GO_terms = ''))
detox.gene.conversion <- gene.table[Gene.name %in% toupper(detox.genes), 
                                    .(Gene.id = Gene_stable_ID, Gene.name = str_to_title(Gene.name))]

# We shrink the modal copy number table to be just the genes we are interested in. These are the detox 
# genes and also some of the ones in the Ace1 deletion in order to get the copy number of that. 
# We also add the cluster Coeaexg
ace1.del.genes <- c('AGAP001358', 'AGAP001360', 'AGAP001361', 'AGAP001362', 'AGAP001363', 'AGAP001364', 'AGAP001365', 'AGAP001366')
genes.of.interest <- c(detox.gene.conversion$Gene.id, ace1.del.genes)
modal.copy.number <- modal.copy.number.all[, c('sample.id', genes.of.interest), with = F] %>%
                     setnames(detox.gene.conversion$Gene.id, detox.gene.conversion$Gene.name) %>%
                     # Add the calculated copy number for Coeaexg
                     merge(phen[good.var.sample.names],
                           by.x = 'sample.id', by.y = 'specimen',
                           all = F
                     ) %>%
                     .[, phenotype := as.factor(phenotype)] %>%
                     setkey(location, insecticide)

source('../shared_functions/R_plotting.r')

# A function that will round to x significant digits, unless the number is small, then just take
# one
signif2 <- function(x, digits, threshold = 0.1){
	small.x <- x < threshold 
	x[small.x] <- signif(x[small.x], 1)
	x[!small.x] <- signif(x[!small.x], digits)
	x
}


target.contable <- function(Dups, include.highvar = T, shrink = T, include.plot = T, 
                            add.sample.size = F, star.small.samples = T, 
                            species.colours = c(coluzzii = 'coral3', gambiae = 'dodgerblue3', arabiensis = 'limegreen'), 
                            text.cell.cex = 0.65, pop.cex = 0.5, dup.cex = 0.5, ...){
	if (include.highvar)
		cnvs.by.population <- population.target.CNVs
	else 
		cnvs.by.population <- population.target.goodvar.CNVs
	population <- droplevels(interaction(cnvs.by.population[, c('location', 'species')]))
	output.table.rownames <- levels(population)
	output.table.colnames <- Dups
	output.table <- matrix(NA, length(output.table.rownames), length(Dups), dimnames = list(output.table.rownames, Dups))
	for (ro in output.table.rownames){
		which.row <- which(population == ro)
		if (length(which.row) == 1)
			output.table[ro, ] <- as.numeric(cnvs.by.population[which.row, Dups])
		else
			output.table[ro, ] <- NA
	}
	if (include.plot){
		if (include.highvar){
			cnvs.by.population <- population.target.CNVs
			popsize <- population.popsize[output.table.rownames]
		}
		else {
			cnvs.by.population <- population.target.goodvar.CNVs
			popsize <- population.goodvar.popsize[output.table.rownames]
		}
		#
		if ((length(names(species.colours)) == 0) | !any(names(species.colours) %in% c('arabiensis', 'gambiae', 'coluzzii', 'intermediate')))
			names(species.colours) <- c('arabiensis', 'gambiae', 'coluzzii', 'intermediate')
		par(mar = c(0,3,3,0), ...)
		plot(c(0,ncol(output.table)), -c(0,nrow(output.table)), type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
		x.matrix <- matrix(1:ncol(output.table), nrow(output.table), ncol(output.table), byrow = T)
		y.matrix <- matrix(1:nrow(output.table), nrow(output.table), ncol(output.table))
		col.matrix <- matrix(rgb(0.8, 0.8, 0.8), nrow(output.table), ncol(output.table))
		text.col.matrix <- c('black', 'white')[(output.table > 0.5) + 1]
		if (add.sample.size){
			text.cell.matrix <- sub('.*\\.', '.', paste(signif2(output.table, 2), ' (', popsize, ')', sep = ''))
			text.cell.matrix[grepl('NA', text.cell.matrix)] <- NA
		}
		else if (star.small.samples){
			text.cell.matrix <- sub('.*\\.', '.', paste(signif2(output.table, 2), c('', ' *', ' **')[(popsize < 5) + (popsize == 1) + 1], sep = ''))
			text.cell.matrix[grepl('NA', text.cell.matrix)] <- NA
		}
		else
			text.cell.matrix <- signif2(output.table, 2)
		for (i in 1:nrow(output.table)){
			this.row <- output.table[i, ]
			if (any(!is.na(this.row))){
				this.species <- sub('^.*\\.', '', rownames(output.table)[i])
				col.matrix[i, !is.na(this.row)] <- lighten.col(species.colours[this.species], sqrt(this.row[!is.na(this.row)]))
			}
		}
		rect(x.matrix - 1, 1 - y.matrix, x.matrix, - y.matrix, col = col.matrix, border = 'white', lwd = text.cell.cex)
		text(x.matrix - 0.5, 0.5 - y.matrix, text.cell.matrix, cex = text.cell.cex, col = text.col.matrix, font = 2)
		text(x.matrix[1,] - 0.5, (par('usr')[4] - par('usr')[3])/100, srt = 45, adj = c(0,0), colnames(output.table), xpd = NA, cex = dup.cex)
		country.names <- sub('\\.[^.]*$', '', rownames(output.table))
		species <- sub('.*\\.', '', rownames(output.table))
		rect(x.matrix[1,1]-1, 1 - y.matrix, par('usr')[1] - (par('usr')[2] - par('usr')[1])/5, -y.matrix, col = sapply(species.colours[species], lighten.col, 0.4), border = NA, xpd = NA)
		text(x.matrix[1,1]-1, 0.5 - y.matrix[,1], paste(country.names, '  \n', species, '  ', sep = ''), adj = 1, cex = pop.cex, xpd = NA)
	}
	if (shrink){
		# We keep the rows where at least one entry is > 0. We then remove the columns where all the remaining
		# rows are NA.
		which.rows <- apply(output.table, 1, sum, na.rm = T) > 0
		output.table <- output.table[which.rows, , drop = F]
		which.columns <- apply(output.table, 2, function(x) all(is.na(x))) == F
		output.table <- output.table[, which.columns, drop = F]
	}
	invisible(output.table)
}

modal.contable <- function(genes, shrink = T, include.plot = T, add.sample.size = F, 
                           star.small.samples = T, use.agaps = F, 
                           species.colours = c(coluzzii = 'coral3', gambiae = 'dodgerblue3', arabiensis = 'limegreen'), 
                           text.cell.cex = 0.65, pop.cex = 0.5, gene.cex = 0.5, ...){
	cnvs.by.population <- population.modal.CNVs
	#
	if (genes[1] %in% gene.table$Gene_stable_ID)
		gene.names <- gene.table[genes, Gene.name]
	else {
		gene.names <- character()
		for (i in 1:length(genes)){
			gene <- genes[i]
			if (toupper(gene) %in% toupper(gene.table$Gene.name)){
				gene.name <- gene.table$Gene.name[toupper(gene.table$Gene.name) == toupper(gene)]
				if (length(gene.name) != 1) 
					stop('There should be only one gene matching the gene name.')
				gene.names[i] <- str_to_title(gene.name)
				genes[i] <- gene.table[Gene.name == gene.name, Gene_stable_ID]
			}
			else
				stop('Could not find the requested gene in gene.table.')
		}
	}
	population <- droplevels(interaction(cnvs.by.population[, c('location', 'species')]))
	output.table.rownames <- levels(population)
	output.table.colnames <- if(use.agaps) genes else gene.names
	output.table <- matrix(NA, length(output.table.rownames), length(output.table.colnames), dimnames = list(output.table.rownames, output.table.colnames))
	gene.has.cnv <- genes %in% colnames(cnvs.by.population)
	for (ro in output.table.rownames){
		which.row <- which(population == ro)
		if (length(which.row) != 1)
			output.table[ro, ] <- NA
		else {
			for (i in 1:ncol(output.table)){
				co <- output.table.colnames[i]
				if (gene.has.cnv[i])
					output.table[ro, co] <- as.numeric(cnvs.by.population[which.row, genes[i]])
				else
					output.table[ro, co] <- 0
			}
		}
	}
	if (include.plot){
		if (is.null(names(species.colours)) | !any(names(species.colours) %in% c('gambiae', 'coluzzii')))
			names(species.colours) <- c('gambiae', 'coluzzii')
		par(mar = c(0,3,3,0), ...)
		plot(c(0,ncol(output.table)), -c(0,nrow(output.table)), type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
		x.matrix <- matrix(1:ncol(output.table), nrow(output.table), ncol(output.table), byrow = T)
		y.matrix <- matrix(1:nrow(output.table), nrow(output.table), ncol(output.table))
		col.matrix <- matrix(rgb(0.8, 0.8, 0.8), nrow(output.table), ncol(output.table))
		popsize <- population.goodvar.popsize[output.table.rownames]
		text.col.matrix <- c('black', 'white')[(output.table > 0.5) + 1]
		if (add.sample.size){
			text.cell.matrix <- sub('.*\\.', '.', paste(signif2(output.table, 2), ' (', popsize, ')', sep = ''))
			text.cell.matrix[grepl('NA', text.cell.matrix)] <- NA
		}
		else if (star.small.samples){
			text.cell.matrix <- sub('.*\\.', '.', paste(signif2(output.table, 2), c('', ' *', ' **')[(popsize < 5) + (popsize == 1) + 1], sep = ''))
			text.cell.matrix[grepl('NA', text.cell.matrix)] <- NA
		}
		else {
			text.cell.matrix <- signif2(output.table, 2)
		}
		for (i in 1:nrow(output.table)){
			this.row <- output.table[i, ]
			if (any(!is.na(this.row))){
				this.species <- sub('^.*\\.', '', rownames(output.table)[i])
				col.matrix[i, !is.na(this.row)] <- lighten.col(species.colours[this.species], sqrt(this.row[!is.na(this.row)]))
			}
		}
		rect(x.matrix - 1, 1 - y.matrix, x.matrix, - y.matrix, col = col.matrix, border = 'white', lwd = text.cell.cex)
		text(x.matrix - 0.5, 0.5 - y.matrix, text.cell.matrix, cex = text.cell.cex, col = text.col.matrix, font = 2)
		text(x.matrix[1,] - 0.5, (par('usr')[4] - par('usr')[3])/100, srt = 45, adj = c(0,0), colnames(output.table), xpd = NA, cex = gene.cex)
		country.names <- sub('\\.[^.]*$', '', rownames(output.table))
		species <- sub('.*\\.', '', rownames(output.table))
		rect(x.matrix[1,1]-1, 1 - y.matrix, par('usr')[1] - (par('usr')[2] - par('usr')[1])/5, -y.matrix, col = sapply(species.colours[species], lighten.col, 0.4), border = NA, xpd = NA)
		text(x.matrix[1,1]-1, 0.5 - y.matrix[,1], paste(country.names, '  \n', species, '  ', sep = ''), adj = 1, cex = pop.cex, xpd = NA)
	}
	if (shrink){
		# We keep the rows where at least one entry is > 0. We then remove the columns where all the remaining
		# rows are NA.
		which.rows <- apply(output.table, 1, sum, na.rm = T) > 0
		output.table <- output.table[which.rows, , drop = F]
		which.columns <- apply(output.table, 2, function(x) all(is.na(x))) == F
		output.table <- output.table[, which.columns, drop = F]
	}
	invisible(output.table)
}


contable <- function(x, ...){
	if (grepl('_Dup', x[1]) | grepl('_Del', x[1]) | grepl('_Any', x[1]))
		target.contable(x, ...)
	else
		modal.contable(x, ...)
}

png('detox_gene_modal_CNVs.png', width = 4.5, height = 0.7, units = 'in', res = 300)
par(family = 'Arial')
contable(detox.genes, 
         text.cell.cex = 0.35,
         pop.cex = 0.35,
         gene.cex = 0.35,
         mai = c(0,0.15,0.2,0)
)
dev.off()

# And a pdf version 
pdf('detox_gene_modal_CNVs.pdf', width = 4.5, height = 0.7)
contable(detox.genes, 
         text.cell.cex = 0.35,
         pop.cex = 0.35,
         gene.cex = 0.35,
         mai = c(0,0.15,0.2,0)
)
dev.off()

png('Ace1_diagnostic_read_CNVs.png', width = 2, height = 0.7, units = 'in', res = 300)
par(family = 'Arial')
contable(Dup.clusters$Ace1, 
         text.cell.cex = 0.35,
         pop.cex = 0.35,
         dup.cex = 0.35,
         mai = c(0,0.25,0.23,0.05)
)
dev.off()

png('Cyp6aap_diagnostic_read_CNVs.png', width = 6.5, height = 0.8, units = 'in', res = 300)
par(family = 'Arial')
contable(Dup.clusters$Cyp6aap, 
         text.cell.cex = 0.35,
         pop.cex = 0.35,
         dup.cex = 0.35,
         mai = c(0,0.09,0.33,0.05)
)
dev.off()

png('Cyp6mz_diagnostic_read_CNVs.png', width = 2.5, height = 0.8, units = 'in', res = 300)
par(family = 'Arial')
contable(Dup.clusters$Cyp6mz, 
         text.cell.cex = 0.35,
         pop.cex = 0.35,
         dup.cex = 0.35,
         mai = c(0,0.23,0.35,0.18)
)
dev.off()

png('Gste_diagnostic_read_CNVs.png', width = 3, height = 0.8, units = 'in', res = 300)
par(family = 'Arial')
contable(Dup.clusters$Gstue, 
         text.cell.cex = 0.35,
         pop.cex = 0.35,
         dup.cex = 0.35,
         mai = c(0,0.21,0.3,0.18)
)
dev.off()

png('Cyp9k1_diagnostic_read_CNVs.png', width = 6, height = 0.8, units = 'in', res = 300)
par(family = 'Arial')
contable(Dup.clusters$Cyp9k1, 
         text.cell.cex = 0.35,
         pop.cex = 0.35,
         dup.cex = 0.35,
         mai = c(0,0.09,0.32,0)
)
dev.off()

png('Coeae_diagnostic_read_CNVs.png', width = 5, height = 0.8, units = 'in', res = 300)
par(family = 'Arial')
contable(c(Dup.clusters$Coeaexf, Dup.clusters$Coeaexg), 
         text.cell.cex = 0.35,
         pop.cex = 0.35,
         dup.cex = 0.35,
         mai = c(0,0.25,0.33,0.18)
)
dev.off()

# We can also plot a histogram of model copy numbers for coeaexg
png('Coeaexg_copy_number_histogram.png', width = 3.5, height = 4, units = 'in', res = 300)
par(mar = c(1.5,2.6,0.5,0.5), mgp = c(1.6,0.4,0), tcl = -0.3, lwd = 1.5, cex = 0.7)
num.bins <- max(modal.copy.number$Coeaexg + 1)
bin.pos <- c(0, seq(5, num.bins, 5))
hist(modal.copy.number$Coeaexg, right = F, main = '', 
     breaks = num.bins, xaxt = 'n', col = lighten.col('#E7B800', alpha = 0.25) , border = '#E7B800')
axis(1, at = bin.pos + 0.5, labels = bin.pos, lwd = 0, mgp = c(0,-0.6,0))
mtext('Copy number', 1, line = 0.5, cex = 0.7)
dev.off()

glm.up <- function(input.table, list.of.markers = markers, rescolumn = 'phenotype', control.for = character(), glm.function = NULL, verbose = T){
	# Check whether the markers and random effects are present in the data.frame
	if (sum(list.of.markers %in% colnames(input.table)) != length(list.of.markers))
		stop('Some of the requested markers were not found in genotypes table.')
	if (sum(control.for %in% colnames(input.table)) != length(control.for))
		stop('At least one random effect was not found in genotypes table.')
	if (!(rescolumn %in% colnames(input.table)))
		stop('Resistance column not found in genotypes table.')
	# Remove any requested random effects that have only 1 level
	random.effects <- character()
	for (this.control in control.for){
		level.counts <- tapply(input.table[,this.control], input.table[,this.control], length)
		if (max(level.counts, na.rm = T) < sum(level.counts, na.rm = T)) 
			random.effects <- c(random.effects, this.control)
		else
			cat('Removing random effect ', this.control, ' was removed because it is invariable.')
	}
	# If you did not set a glm.function, decide which glm function you are going to use, based on whether mixed 
	# modelling will be necessary
	if (is.null(glm.function))
		glm.function <- ifelse(length(random.effects) > 0, 'glmer', 'glm')
	if (glm.function == 'glm'){
		# Create the random effect string, which will be empty if we have no random effects
		random.effects.string <- ifelse(length(random.effects) > 0, paste(' +', paste(random.effects, collapse = ' + ', sep = '')), '')
		# Set the name of the column containing the P.value in the anova function (which is different depending
		# on the glm function you use
		P.val.column <- 'Pr(>Chi)'
	}
	else if (glm.function %in% c('glmer', 'glmmTMB')){
		random.effects.string <- ifelse(length(random.effects) > 0, paste(' +', paste('(1|', random.effects, ')', collapse = ' + ', sep = '')), '')
		P.val.column <- 'Pr(>Chisq)'
	}
	if (verbose)
		cat('Using the following string to control for confounding factors: ', random.effects.string, '\n', sep = '')
	# We remove markers for which there is no variation in the dataset or for which some alleles are too rare. 
	if (verbose)
		cat('\nDetecting invariable and nearly invariable markers.\n')
	kept.markers <- character()
	invariable.markers <- character()
	for (this.marker in list.of.markers){
		allele.counts <- tapply(input.table[,this.marker], input.table[,this.marker], length)
		if (max(allele.counts, na.rm = T) <= (sum(allele.counts, na.rm = T) - 2)) 
			kept.markers <- c(kept.markers, this.marker)
		else
			invariable.markers <- c(invariable.markers, this.marker)
	}
	if (length(kept.markers) == 0){
		cat('\nFail. None of the markers provided were sufficiently variable.\n')
		return(NULL)
	}
	# We check whether there are any ordered factors and recode them as numeric
	converted.table <- input.table
	has.ordered.factors <- F
	for (this.marker in kept.markers){
		if ('ordered' %in% class(converted.table[, this.marker])){
			if (verbose)
				cat('Converting ordered factor ', this.marker, ' to numeric.\n', sep = '')
			converted.table[, this.marker] <- as.numeric(converted.table[, this.marker])
			has.ordered.factors <- T
		}
	}
	# We do the glm analysis directly on the table from the global environment rather than the argument, this 
	# way the table that was used is recorded in the output. If we had to convert the ordered factors, then
	# we are forced to use a new table
	if (has.ordered.factors){
		working.table.name <- make.names(paste(deparse(substitute(input.table)), '_numeric_conversion', sep = ''))
		eval(parse(text = paste(working.table.name, '<- converted.table')))
	}
	else
		working.table.name <- deparse(substitute(input.table))
	
	# For each marker, we calculate its pseudo-R2 and P-value compared to the null model.
	if (verbose)
		cat('\nAnalysing markers independently.\n')
	individual.markers <- data.frame(P = numeric(), pseudo.R2 = numeric())
	for (this.marker in kept.markers){
		# Remove the Na values for this marker
		this.table <- converted.table[!is.na(converted.table[,this.marker]),]
		# Build the model 
		this.model <- eval(parse(text = paste(glm.function, '(', rescolumn, ' ~ ', this.marker, random.effects.string, ', data = this.table, family = binomial)', sep = '')))
		# Build the null model
		this.null.model <- eval(parse(text = paste(glm.function, '(', rescolumn, ' ~ 1 ', random.effects.string, ', data = this.table, family = binomial)', sep = '')))
		# Get the stats
		this.p <- anova(this.model, this.null.model, test = 'Chisq')[[P.val.column]][2]
		# Report pseudo Rsquared if we used GLM and if we have the modEvA package
		if (('modEvA' %in% (.packages())) & (glm.function == 'glm'))
			this.pseudo.r2 <- mean(unlist(RsqGLM(this.model)))
		else
			this.pseudo.r2 <- NA
		individual.markers[this.marker, ] <- c(this.p, this.pseudo.r2)
	}
	
	# We now build the null model and add markers one by one until all markers are significant
	working.markers <- character()
	# We'll keep track of the markers with perfect correlation
	correlated.markers <- character()
	if (verbose)
		cat('\nRunning commentary on model optimisation:\n\n')
	while(length(working.markers) < length(kept.markers)){
		# Build the model using the working.markers
		if (length(working.markers)){
			old.model.text <- paste(glm.function, '(', rescolumn, ' ~ ', paste(working.markers, collapse = ' + '), random.effects.string, ', data = ', working.table.name, ', family = binomial)', sep = '')
			old.model <- eval(parse(text = old.model.text))
			# Check the remaining markers as some of them may become monomorphic when NAs from the current marker are 
			# taken into account
			for (this.marker in setdiff(kept.markers, working.markers)){
				markers.subset.genotypes <- input.table[, c(working.markers, this.marker), drop = F]
				markers.subset.genotypes <- markers.subset.genotypes[complete.cases(markers.subset.genotypes), , drop = F]
				number.of.alleles <- unique(markers.subset.genotypes[ ,this.marker])
				if (length(number.of.alleles) < 2) {
					if (verbose){
						cat('Removing marker ', this.marker, ' as it has no variation left when NAs at previously added ',
							'loci are removed.\n\n', sep = '')
					}
					kept.markers <- setdiff(kept.markers, this.marker)
				}
			}
			# If we have removed all the remaining markers, we quit the loop and report the final model. 
			if (length(working.markers) == length(kept.markers)){
				final.model <- old.model
				if (verbose){
					cat('\tNo further markers are significant, keeping final model:\n')
					print(final.model)
				}
				break
			}
		}
		else{
			old.model.text <- paste(glm.function, '(', rescolumn, ' ~ 1 ', random.effects.string, ', data = ', working.table.name, ', family = binomial)', sep = '')
			old.model <- eval(parse(text = old.model.text))
		}
		if (verbose)
			cat('Building model:', old.model.text, '\n')
		p.values <- numeric()
		for (this.marker in setdiff(kept.markers, working.markers)){
			new.model <- eval(parse(text = paste('update(old.model, .~.+', this.marker, ')', sep = '')))
			# Check that the new model doesn't have fewer rows than the old one (if the added marker had unique NAs)
			if (length(fitted(new.model)) == length(fitted(old.model))){
				this.p.value <- anova(old.model, new.model, test = 'Chisq')[[P.val.column]][2]
			}
			# Otherwise, we need to rebuild the models with those samples removed.
			else{
				if (has.ordered.factors)
					reduced.input.table <- converted.table[!is.na(converted.table[,this.marker]),]
				else 
					reduced.input.table <- input.table[!is.na(input.table[,this.marker]),]
				if (length(working.markers))
					temp.old.model <- eval(parse(text = paste(glm.function, '(', rescolumn, ' ~ ', paste(working.markers, collapse = ' + '), random.effects.string, ', data = reduced.input.table, family = binomial)', sep = '')))
				else 
					temp.old.model <- eval(parse(text = paste(glm.function, '(', rescolumn, ' ~ 1 ', random.effects.string, ', data = reduced.input.table, family = binomial)', sep = '')))
				new.model <- eval(parse(text = paste('update(temp.old.model, .~.+', this.marker, ')', sep = '')))
				this.p.value <- anova(temp.old.model, new.model, test = 'Chisq')[[P.val.column]][2]
			}
			# If the p.value was NA, then there is perfect correlation or something else was wrong. Set the p-value 
			# to Inf and move on to the next marker
			if (is.na(this.p.value)){
				if (verbose)
					cat('\tCould not calculate p-value for marker ', this.marker, '.\n\n', sep = '')
				p.values[this.marker] <- Inf
			}
			else{
				p.values[this.marker] <- this.p.value
			}
		}
		if (min(p.values) <= 0.05){
			# Add the lowest significant p-value
			if (verbose)
				cat('\tAdding marker ', names(p.values)[which.min(p.values)], ' as the lowest significant marker (P = ', min(p.values), ').\n\n', sep = '')
			marker.to.add <- names(p.values)[which.min(p.values)]
			working.markers <- c(working.markers, marker.to.add)
			if (length(working.markers) == length(kept.markers)){
				# If all markers have been added, then we have the final model
				final.model <- new.model
				if (verbose)
					cat('\tNo markers left to add.\n')
			}
		}
		else {
			final.model <- old.model
			if (verbose){
				cat('\tNo further markers are significant, keeping final model:\n')
				print(final.model)
			}
			break
		}
	}
	if (verbose){
		if (length(working.markers) == 0)
			cat('Final model was the null model.\n\n')
		else 
			cat('Final model contained ', length(working.markers), ' parameters: ', paste(working.markers, collapse = ','), '.\n\n', sep = '')
	}
	# Now get the p-values and pseudo R-squared value for all the variables in the final model, when added as the 
	# last variable
	if (length(working.markers) > 0){
		deviance.effect <- numeric()
		final.p.values <- numeric()
		for (this.marker in working.markers){
			# Remove the Na values for this marker
			this.table <- converted.table[!is.na(converted.table[,this.marker]),]
			# Build the model 
			reduced.final.model <- update(final.model, data = this.table)
			reduced.model <- eval(parse(text = paste('update(reduced.final.model, .~.-', this.marker, ')', sep = '')))
			# Get the stats
			dev.eff <- deviance(reduced.model) - deviance(final.model)
			deviance.effect[this.marker] <- ifelse(is.null(dev.eff), NA, dev.eff)
			final.p.values[this.marker] <- anova(reduced.model, reduced.final.model, test = 'Chisq')[[P.val.column]][2]
		}
	}
	else {
		final.p.values <- NA
		deviance.effect <- NA
	}
	cat('\n')
	#
	final.model.sig <- data.frame(P = final.p.values, deviance = deviance.effect)
	print(final.model.sig)
	list('invariable.markers' = invariable.markers, 'correlated.markers' = correlated.markers, 'sig.alone' = individual.markers, 'final.model' = final.model, 'final.sig' = final.model.sig)
}

# Let's test things using copy number in each population independently
cat('\n\t#######################')
cat('\n\t# Copy number testing #')
cat('\n\t#######################\n')

detox.genes.with.cnvs <- intersect(detox.gene.conversion$Gene.id, names(modal.CNV.table))
dtx.gene.conversion <- detox.gene.conversion[Gene.id %in% detox.genes.with.cnvs]
detox.gene.freq <- aggregate(modal.CNV.table[, dtx.gene.conversion$Gene.id, with = F],
                             phen[modal.CNV.table$sample.id, .(location, insecticide)], 
                             function(x) sum(x) / length(x)
                   ) %>%
                   setnames(dtx.gene.conversion$Gene.id, dtx.gene.conversion$Gene.name)
rownames(detox.gene.freq) <- paste(detox.gene.freq$location, detox.gene.freq$insecticide, sep = '.')
detox.gene.freq <- detox.gene.freq[sort(dtx.gene.conversion$Gene.name)]

for (insecticide in c('Delta', 'PM')){
	for (location in c('Moshi', 'Muleba')){
		if (insecticide == 'PM' & location == 'Muleba')
			next
		# For some reason, this doesn't work if I do the indexing directly. I have to create this object first
		index <- list(location, insecticide)
		# Keep genes where CNV exists in at least 5% of samples. 
		genes.to.model <- (detox.gene.freq[paste(index, collapse = '.'), ] >= 0.05) %>%
		                  colnames(.)[.]
		cat('\n\tModal copy number test for ', paste(index, collapse = '_'), ':\n\n', sep = '')
		cat('Analysing genes with CNV at >= 5% freq (', paste(genes.to.model, collapse = ', '), ').\n\n', sep = '')
		test.table <- as.data.frame(modal.copy.number[index])
		print(glm.up(test.table, genes.to.model, 'phenotype'))
	}
}

# Let's try a glm that includes population as a random factor. 
genes.to.model <- c('Cyp6aa1', 'Cyp6z2', 'Cyp6p3','Gste2', 'Cyp9k1', 'Coeae2f', 'Coeaexg')
delta.table <- as.data.frame(modal.copy.number[insecticide == 'Delta'])
cat('\n\nAll populations Delta:\n')
delta.model <- glm.up(delta.table, genes.to.model, 'phenotype', control.for = 'location', glm.function = 'glmmTMB')
print(delta.model)
# Cyp6aa1 is significant

cat('\n\t##############################')
cat('\n\t# Cyp6aap Dup allele testing #')
cat('\n\t##############################\n')

# For Cyp6aap, let's look at each Dup in turn. We'll do that by population because there are different CNVs in
# different populations. We'll use CNVs at a minimum of 5% freq
# Moshi:
cat('\n\nCyp6aap alleles in Moshi:\n')
moshi.delta.samples <- phen[location == 'Moshi', specimen]
moshi.delta.allele.table <- as.data.frame(target.CNV.table[moshi.delta.samples])
moshi.delta.model <- glm.up(moshi.delta.allele.table, c(paste('Cyp6aap_Dup', 31:33, sep = ''), 'Cyp9k1_Dup18'), 'phenotype')
print(moshi.delta.model)
# Dup33 is significant

# Muleba:
cat('\n\nCyp6aap alleles in Muleba:\n')
muleba.delta.samples <- phen[location == 'Muleba' & insecticide == 'Delta', specimen]
muleba.delta.allele.table <- as.data.frame(target.CNV.table[muleba.delta.samples])
muleba.delta.model <- glm.up(muleba.delta.allele.table, c('Cyp6aap_Dup33'), 'phenotype')
print(muleba.delta.model)
# Again, Dup33 is significant

# Now do that again, but include the presence / absence of the sweep:
cat('\n\nCyp6aap alleles in Moshi, controlling for sweep:\n')
moshi.haps <- fread('../haplotypes/Moshi_arabiensis_Delta_Cyp6.csv', key = 'sample_name')
moshi.delta.allele.table <- cbind(moshi.delta.allele.table, moshi.haps[moshi.delta.allele.table$sample.id, -c(1,2)])
moshi.delta.haps.model <- glm.up(moshi.delta.allele.table, c('Cyp6aap_Dup33', 'cluster_1'), 'phenotype')
# Still just Dup33

cat('\n\nCyp6aap alleles in Muleba, controlling for sweep:\n')
muleba.haps <- fread('../haplotypes/Muleba_arabiensis_Delta_Cyp6.csv', key = 'sample_name')
muleba.delta.allele.table <- cbind(muleba.delta.allele.table, muleba.haps[muleba.delta.allele.table$sample.id, -c(1,2)])
muleba.delta.haps.model <- glm.up(muleba.delta.allele.table, c('Cyp6aap_Dup33', 'cluster_1', 'cluster_2'), 'phenotype')
# Still just Dup33

# Write tables to file
# Find the Dups found at least once in the data
Dups.present <- known.Dups[apply(target.CNV.table[, ..known.Dups], 2, sum) > 0]
output.table <- merge(phen[, .(specimen, species, location, insecticide)],
                      target.CNV.table[, c('sample.id', 'phenotype', 'High.var.sample', ..Dups.present)], 
                      by.x = 'specimen',
                      by.y = 'sample.id',
                      all.x = F,
                      all.y = T) %>%
                setnames(., 'specimen', 'sample.id') %>%
                merge(., modal.copy.number[, c('sample.id', detox.genes), with = F], by = 'sample.id', all = T) %>%
				.[order(location, insecticide)]

fwrite(data.table(sample_id = names(coeaexg.median.cn), copy_number = coeaexg.median.cn), 'Coeaexg_copy_number.csv', sep = '\t')
fwrite(output.table, 'CNV_calls.csv', sep = '\t')


