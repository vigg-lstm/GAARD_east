library(data.table)

meta <- fread('../data/tanzania.samples.meta.csv', key = 'sample_id')
meta.sample.names <- meta$partner_sample_id
phen <- fread('../data/sample_phenotypes_EA.csv', key = 'specimen')[sort(meta.sample.names), ]

# Load the sib groups information
sib.groups <- fread('../NGSrelate/full_relatedness_tanzania/sib_group_table.csv')
sibs.to.remove <- sib.groups[keep == F, sample.name]
males <- meta[sex_call == 'M', partner_sample_id]
samples.to.exclude = readLines('../NGSrelate/full_relatedness_tanzania/samples_to_exclude.csv')
samples.to.remove <- c(sibs.to.remove, males, samples.to.exclude)

sample.sizes <- phen[, .(dead = sum(phenotype == 'dead'), alive = sum(phenotype == 'alive')), 
				       keyby = c('location', 'insecticide')]

nomales.phen <- phen[!(males)]
nomales.sample.sizes <- nomales.phen[, .(dead = sum(phenotype == 'dead'), alive = sum(phenotype == 'alive')), 
				                       keyby = c('location', 'insecticide')]

final.phen <- phen[!(samples.to.remove)]
final.sample.sizes <- final.phen[, .(dead = sum(phenotype == 'dead'), alive = sum(phenotype == 'alive')), 
				                   keyby = c('location', 'insecticide')]

cat('\nSample sizes before removing males and sibs:\n\n')
print(sample.sizes)

cat('\nSample sizes after removing males:\n\n')
print(nomales.sample.sizes)

cat('\n\nSample sizes after removing males and sibs:\n\n')
print(final.sample.sizes)

fwrite(nomales.sample.sizes, file = 'sample_sizes_tanzania.csv', sep = '\t')
fwrite(final.sample.sizes, file = 'final_sample_sizes_tanzania.csv', sep = '\t')

