library(data.table)
library(stringr)

arg.values <- commandArgs(trailingOnly=T)

meta.fn <- arg.values[1]
output.fn <- arg.values[2]
{if (length(arg.values) == 2)
	samples.to.remove <- 'none'
else
	samples.to.remove <- readLines(arg.values[3])
}

arg.values <- commandArgs(trailingOnly=T)
cat('Running with metadata:', meta.fn, 
    'and output filename:', output.fn, 
    'samples to remove:', paste(samples.to.remove, sep = ','), '\n', sep = '\n')

# Load the phenotypes (original GAARD metadata) and the VObs metadata (provided by Sanger, contains only 
# sequenced samples. 
phenotypes <- fread('../data/sample_phenotypes_EA.csv', key = 'specimen')
vobs.meta <- fread(meta.fn, key = 'partner_sample_id')
# Remove samples that need excluding
samples.to.include <- setdiff(
	vobs.meta[sex_call == 'F', partner_sample_id],
	samples.to.remove
)

# Reduced the phenotype table to only include sequenced samples
phenotypes <- phenotypes[samples.to.include]

# Have a column to indicate population (location x species x insecticide)
phenotypes$population <- with(phenotypes, paste(location, species, insecticide, sep = '.'))

# Reorder and remove columns
column.order <- c('specimen', 'country', 'location', 'species', 'insecticide', 'population', 'phenotype')
phenotypes <- phenotypes[, ..column.order]

# A function that will return a list of length k, where each k is a random shuffling of the input vector
shuffle <- function(x, k)
	data.table(replicate(k, sample(x, length(x), replace = F)))

# Create 1000 randomisations of the phenotype labels, stratified by population, and add them to the phenotypes
# table
set.seed(42)
num.randomisations <- 10000
replicate.names <- paste('r', str_pad(1:num.randomisations, nchar(as.character(num.randomisations)), pad = 0), sep = '')
phenotypes[, (replicate.names) := shuffle(phenotype, num.randomisations), by = population]

# Write the table to file
fwrite(phenotypes, output.fn, sep = '\t', quote = F)

