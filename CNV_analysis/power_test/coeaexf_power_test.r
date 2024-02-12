library(data.table)
library(magrittr)
library(stringr)
library(glmmTMB)
library(abind)

study.ids.table <- fread('../../data/study_ids.csv')
study.ids.table[, study_id := paste('v3.7_', study_id, sep = '')]
meta <- fread('../../data/tanzania.samples.meta.csv', key = 'sample_id')
meta.sample.names <- meta$partner_sample_id
phen <- fread('../../data/sample_phenotypes_EA.csv', key = 'specimen')[sort(meta.sample.names), ]

# Identify samples for removal (sibs, males and samples to exlude due to unusual
# relatedness)
sib.groups <- fread('../../NGSrelate/full_relatedness_tanzania/sib_group_table.csv')
sibs.to.remove <- sib.groups[keep == F, sample.name]
males <- meta[sex_call == 'M', partner_sample_id]
samples.to.exclude = readLines('../../NGSrelate/full_relatedness_tanzania/samples_to_exclude.csv')
samples.to.remove <- c(sibs.to.remove, males, samples.to.exclude)

# Now get the modal CNVs by gene
load.modal.cnv.table <- function(study.id){
	modal.cnv.table <- fread(paste('../Ag1000G_CNV_data', study.id, 'modal_CNVs/modal_copy_number_arabiensis.csv', sep = '/'))
	these.sample.names <- meta[modal.cnv.table$V1, partner_sample_id]
	modal.cnv.table$V1 <- these.sample.names
	colnames(modal.cnv.table)[1] <- 'sample.id'
	modal.cnv.table
}

modal.copy.number.all <- lapply(unique(study.ids.table$study_id), load.modal.cnv.table) %>%
                         rbindlist() %>%
                         .[!(sample.id %in% samples.to.remove)] %>%
                         .[high_variance == F, -c('sex_call', 'high_variance')]

detox.genes <- 'Coeae1f'

# Get the names of the detox genes
gene.table <- fread('../Ag1000G_CNV_data/gene_annotation_fullgenetable.csv', key = 'Gene_stable_ID', check.names = T)
# Coeae1f is absent from the gff, so add it manually
gene.table['AGAP006227', Gene.name := 'COEAE1F']
# We also add a special entry for the Coeaexg cluster, for subsequent functions to work
detox.gene.conversion <- gene.table[Gene.name %in% toupper(detox.genes), 
                                    .(Gene.id = Gene_stable_ID, Gene.name = str_to_title(Gene.name))]

# We shrink the modal copy number table to be just the genes we are interested in. These are the detox 
# genes and also some of the ones in the Ace1 deletion in order to get the copy number of that. 
# We also add the cluster Coeaexg
genes.of.interest <- detox.gene.conversion$Gene.id
modal.copy.number <- modal.copy.number.all[, c('sample.id', genes.of.interest), with = F] %>%
                     setnames(detox.gene.conversion$Gene.id, detox.gene.conversion$Gene.name) %>%
                     # In the original version of the script, we indexed phen by good.var.sample.names, but we don't
                     # actually need to do this, because modal.copy.number.all has already been filtered, and we have
                     # set all = F in the merge function
                     merge(phen,
                           by.x = 'sample.id', by.y = 'specimen',
                           all = F
                     ) %>%
                     .[, phenotype := as.factor(phenotype)] %>%
                     setkey(location, insecticide)

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
	if (verbose)
		cat('\n')
	#
	final.model.sig <- data.frame(P = final.p.values, deviance = deviance.effect)
	if (verbose)
		print(final.model.sig)
	list('invariable.markers' = invariable.markers, 'correlated.markers' = correlated.markers, 'sig.alone' = individual.markers, 'final.model' = final.model, 'final.sig' = final.model.sig)
}

set.seed(42)
index <- list('Moshi', 'PM')
model.table <- as.data.frame(modal.copy.number[index, c('Coeae1f', 'phenotype')])

# Here is the contingency table for Coeae1f from Obuasi, converted to presence / absence (only one sample in 
# Obuasi has copy number > 1)
wa.equivalent.table <- matrix(c(43,10,34,2), 2, 2)
# Calculate the proportion of alive individuals in each of the wt and CNV-carrying sub-groups
ps <- wa.equivalent.table[,1] / apply(wa.equivalent.table, 1, sum)
p1 <- ps[1]
p2 <- ps[2]

# Use this to simulate data with the number of wt and CNV-carrying samples in Moshi
k = 1000
num.wt <- sum(model.table$Coeae1f == 0)
num.mut <- sum(model.table$Coeae1f > 0)
results <- numeric(k)
for (i in 1:k){
	simulated.wt <- sample(c('dead','alive'), num.wt, replace = T, prob = c(1-p1, p1))
	simulated.mut <- sample(c('dead','alive'), num.mut, replace = T, prob = c(1-p2, p2))
	simulated.table <- rbind(data.frame(Coeae1f = rep(0, num.wt), phenotype = factor(simulated.wt, levels = c('alive', 'dead'))),
							 data.frame(Coeae1f = rep(1, num.mut), phenotype = factor(simulated.mut, levels = c('alive', 'dead')))
					   )
	results[i] <- glm.up(simulated.table, 'Coeae1f', 'phenotype', verbose = F)$sig.alone$P
}


# What was our statistical power?
stat.power <- sum(results < 0.05) / k
cat('The statistical power was', stat.power, '\n')

