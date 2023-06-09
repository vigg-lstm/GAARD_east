library(data.table)
library(magrittr)
library(stringr)

# Load meta data
meta <- fread('Tanzania_metadata.csv')
phenotypes <- fread('../../data/sample_phenotypes_EA.csv')
meta <- merge(meta, phenotypes[, .(specimen, insecticide, phenotype)], 
              by.x = 'partner_sample_id', by.y = 'specimen',
			  all.x = T, all.y = F
)

# Load relatedness results
ngs.relate.output.fn <- 'tanzania.allsnps.king.csv'
ngs.relate <- fread(ngs.relate.output.fn, sep = '\t')
# Change the name of the king column for easier typing
colnames(ngs.relate)[3] <- tolower(colnames(ngs.relate)[3])

# Change numerical indices to sample names. Indices are 0-indexed, so add +1
ngs.relate$a <- meta[ngs.relate$a + 1, 'partner_sample_id']
ngs.relate$b <- meta[ngs.relate$b + 1, 'partner_sample_id']

# Now set the metadata table key to be the GAARD sample name. 
setkey(meta, 'partner_sample_id')

ngs.relate$species.a <- meta[ngs.relate$a, aim_species_gambcolu_arabiensis]
ngs.relate$species.b <- meta[ngs.relate$b, aim_species_gambcolu_arabiensis]
ngs.relate$location.a <- meta[ngs.relate$a, location]
ngs.relate$location.b <- meta[ngs.relate$b, location]
ngs.relate$insecticide.a <- meta[ngs.relate$a, insecticide]
ngs.relate$insecticide.b <- meta[ngs.relate$b, insecticide]
ngs.relate$phenotype.a <- meta[ngs.relate$a, phenotype]
ngs.relate$phenotype.b <- meta[ngs.relate$b, phenotype]

# Overall histogram 
png('king_histogram.png', width = 6, height = 2.2, units = 'in', res = 300)
par(cex = 0.4, lwd = 0.6, mar = c(3,2,2,0), oma = c(0,1,0,0), mgp = c(2,0.8,0), xpd = NA)
hist(ngs.relate$king, col = rgb(0,0,0.8), breaks = 100, xlab = 'King', ylab = '')
# Within the histogram, show the proportion of each bar that are from same location 
with(ngs.relate[location.a == location.b], hist(king, col = rgb(0.6,0.6,1), add = T, breaks = 100))
dev.off()

# Mode of relatedness is above 0. 

# Now zoom in to try and see sib clusters. 
png('king_histogram_zoomed.png', width = 4, height = 2.2, units = 'in', res = 300)
par(cex = 0.4, lwd = 0.6, mar = c(3,2,2,0), oma = c(0,1,0,0), mgp = c(2,0.8,0))
hist(ngs.relate$king, col = rgb(0,0,0.8), breaks = 200, xlab = 'King', ylab = '', ylim = c(0,100))
mtext('frequency', 2, 0, T, cex = 0.45)
with(ngs.relate[location.a == location.b], hist(king, col = rgb(0.6,0.6,1), add = T, breaks = 200))
dev.off()
# There are some extra modes there, very low, which include cross-location pairs. But the really 
# highly-related pairs (possible sibs) are all same-pop. 

# From playing around with the data, sample 'EA18-0330' came out as being highly related to many samples. 
# Not sure why this would be or how it could happen. Let's highlight the pairs involving that sample
#png('king_histogram_0330.png', width = 4, height = 2.2, units = 'in', res = 300)
par(cex = 0.4, lwd = 0.6, mar = c(3,2,2,0), oma = c(0,1,0,0), mgp = c(2,0.8,0))
hist1 <- hist(ngs.relate$king, col = rgb(0,0,0.8), breaks = 200, xlab = 'King', ylab = '', ylim = c(0,100))
mtext('frequency', 2, 0, T, cex = 0.45)
with(ngs.relate[(a == 'EA18-0330') | (b == 'EA18-0330')], hist(king, col = rgb(1,0.6,0.6), add = T, breaks = hist1$breaks))
#dev.off()

# So that sample is just really highly related to everything somehow. I wonder if those other modes
# are similarly driven by a single sample. Answer is yes:
odd.samples <- c('EA18-0330', 'EA18-0268', 'EA18-0195', 'EA18-0234')
png('king_histogram_oddsamples.png', width = 4, height = 2.2, units = 'in', res = 300)
par(cex = 0.4, lwd = 0.6, mar = c(3,2,2,0), oma = c(0,1,0,0), mgp = c(2,0.8,0))
hist1 <- hist(ngs.relate$king, col = rgb(0,0,0.8), breaks = 200, xlab = 'King', ylab = '', ylim = c(0,100))
mtext('frequency', 2, 0, T, cex = 0.45)
with(ngs.relate[(a %in% odd.samples) | (b %in% odd.samples)], hist(king, col = rgb(1,0.6,0.6), add = T, breaks = hist1$breaks))
with(ngs.relate[(a %in% odd.samples[2:4]) | (b %in% odd.samples[2:4])], hist(king, col = rgb(0.6,1,0.6), add = T, breaks = hist1$breaks))
with(ngs.relate[(a %in% odd.samples[3:4]) | (b %in% odd.samples[3:4])], hist(king, col = rgb(1,1,0.6), add = T, breaks = hist1$breaks))
with(ngs.relate[(a %in% odd.samples[4]) | (b %in% odd.samples[4])], hist(king, col = rgb(1,0.6,1), add = T, breaks = hist1$breaks))
dev.off()

# Those four samples together pretty much explain all the weirdly high relatedness values. So let's 
# just chuck them out of the dataset. 

ngs.relate.clean <- ngs.relate[!((a %in% odd.samples) | (b %in% odd.samples))]

# A function to build sib groups based on a given kinship threshold value. 
find.sib.groups <- function(sib.table, verbose = F){
	# First, create a list where each entry is a sib pair declared by the sib.table
	pair.list <- lapply(split(sib.table[, .(a, b)], 1:nrow(sib.table)), unlist)
	# Now merge any two entries in that list that share a sample
	for (i in length(pair.list):2){
		shared.samples <- which(sapply(pair.list[1:(i-1)], function(x) any(pair.list[[i]] %in% x)))
		if (length(shared.samples) >= 1){
			pair.list[[shared.samples[1]]] <- unique(c(pair.list[[shared.samples[1]]], pair.list[[i]]))
			pair.list[[i]] <- NULL
		}
	}
	sib.groups <- lapply(pair.list, function(group) sib.table[a %in% group | b %in% group, ])
	cat(length(sib.groups), 'sib groups were found.\n')
	
	cross.study.sibs.per.group <- sapply(sib.groups, function(sib) sum(sib$location.a != sib$location.b))
	population.bridge <- cross.study.sibs.per.group > 0
	if (any(population.bridge))
		cat('\t', sum(population.bridge), ' sib groups bridged across populations.')
	total.cross.study.pairs <- sum(cross.study.sibs.per.group)
	
	# Now go through and deal with groups that don't have coherent sets of sibs
	missing.pair.tables <- list()
	total.missing.pairs <- 0
	for (sib in sib.groups){
		members <- with(sib, unique(unlist(c(a, b))))
		expected.pairs <- combn(members, 2)
		num.missing <- ncol(expected.pairs) - nrow(sib)
		if (num.missing > 0){
			total.missing.pairs <- total.missing.pairs + num.missing
			# Get the table of missing pairs from the full table
			missing.pairs <- c()
			for (j in 1:ncol(expected.pairs)){
				this.pair <- expected.pairs[, j]
				pair.rows <- with(sib, (a %in% this.pair) & (b %in% this.pair))
				# If none of those values are true, then this is a missing pair and needs to be found in the full table
				if (!any(pair.rows))
					missing.pairs <- c(missing.pairs, which(with(ngs.relate, (a %in% this.pair) & (b %in% this.pair))))
			}
			missing.pair.table <- ngs.relate[missing.pairs, ]
			if (verbose)
				cat('\t', num.missing, ' missing pairs were found in group containing ', 
					paste(members, collapse = ', '), 
					' with min, mean and max KING scores of ', 
					paste(round(summary(missing.pair.table$king)[c('Min.', 'Mean', 'Max.')], 2), collapse  =', '),
					' respectively.\n', sep = '')
			missing.pair.tables <- c(missing.pair.tables, list(missing.pair.table))
		}
		else {
			if (verbose)
				cat('\tThere were no missing pairs.\n')
		}
	}
	total.pairs <- sum(sapply(sib.groups, nrow))
	list(sib.groups = sib.groups, 
	     total.pairs = total.pairs, 
	     total.missing.pairs = total.missing.pairs, 
		 missing.pair.ratio = total.missing.pairs / total.pairs,
		 incorrect.group.ratio = length(missing.pair.tables) / length(sib.groups),
	     total.cross.study.pairs = total.cross.study.pairs, 
	     missing.pair.tables = missing.pair.tables, 
	     cross.study.sibs.per.group = cross.study.sibs.per.group)
}

test.threshold <- function(king.thresh, pairs.table, verbose = F){
	cat('\nRunning find.sib.groups with a king threshold of', king.thresh, '\n')
	sib.pairs <- pairs.table[king >= king.thresh, ]
	invisible(find.sib.groups(sib.pairs, verbose = verbose))
}

# Find the "best" threshold to define siblings. 
thresholds <- seq(0.15, 0.3, 0.005)
threshold.exploration <- lapply(thresholds, test.threshold, rbind(ngs.relate.clean))
num.cross.pairs <- sapply(threshold.exploration, function(x) x$total.cross.study.pairs)
total.missing.pairs <- sapply(threshold.exploration, function(x) x$total.missing.pairs)
total.pairs <- sapply(threshold.exploration, function(x) x$total.pairs)
missing.pair.ratio <- sapply(threshold.exploration, function(x) x$missing.pair.ratio)
incorrect.group.ratio <- sapply(threshold.exploration, function(x) x$incorrect.group.ratio)

x11()
par(mfrow = c(1,2))
plot(thresholds, total.pairs, pch = 21, col = 'grey50', bg = 'black', ylim = c(0, max(c(total.pairs, total.missing.pairs))), main = 'same species pairs')
points(thresholds, total.missing.pairs, pch = 21, col = 'grey50', bg = 'red')
points(thresholds, num.cross.pairs, pch = 21, col = 'grey50', bg = 'magenta')
plot(thresholds, missing.pair.ratio, pch = 21, col = 'grey50', bg = rgb(1,0.5,0.5, 0.5), ylim = c(0, max(c(missing.pair.ratio, incorrect.group.ratio))))
points(thresholds, incorrect.group.ratio, col = 'grey50', pch = 21, bg = rgb(0.5,1,0.5, 0.5))
thresholds[order(missing.pair.ratio)]

# So we settle on a threshold of 0.185. This is the highest value before we end up with a missing pair. 
# But actually, in this case there is no difference between using this threshold and the king 
# recommended value of 0.175
best.threshold.results <- threshold.exploration[[which(thresholds == 0.185)]]

# Create a list of sib groups, which includes the experiment that each individual was taken from
create.sib.group.table <- function(sib.pair.table, sib.group.id){
	sib.group.table <- data.table(sample.name = c(sib.pair.table$a, sib.pair.table$b), 
	                              location = c(sib.pair.table$location.a, sib.pair.table$location.b),
	                              species = c(sib.pair.table$species.a, sib.pair.table$species.b),
	                              insecticide = c(sib.pair.table$insecticide.a, sib.pair.table$insecticide.b),
	                              sib.group.id = sib.group.id
	                             ) %>%
	                   {.[!duplicated(.)]} %>%
	                   {.[order(sample.name)]}
	sib.group.table
}

sib.group.table <- mapply(create.sib.group.table, 
                          best.threshold.results$sib.groups, 
                          1:length(best.threshold.results$sib.groups),
                          SIMPLIFY = F) %>%
                   rbindlist()

set.seed(42)
sib.group.table[, keep := sample(c(T, rep(F, nrow(.SD) -1)), nrow(.SD)), by = c('location', 'insecticide', 'sib.group.id')]

cat('\nWe have removed', table(sib.group.table$keep)['FALSE'], 'samples as sibs, broken down as follows:\n\n')
meta[sib.group.table[!(keep), sample.name]] %>%
with(., table(location, phenotype, insecticide)) %>%
print

fwrite(sib.group.table, file = 'sib_group_table.csv', sep = '\t')

