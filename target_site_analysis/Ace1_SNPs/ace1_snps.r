library(data.table)
library(stringr)

load('~/data/ML/GAARD_SNP/Moshi_arabiensis_PM/classical_analysis_sibgroups/classical_analysis_nosibs.Rdata')

ace1.snps <- snp.table[Chrom == '2R' & Pos > 3484100 & Pos < 3495600]

ace1.allele.counts.filename <- 'ace1_allele_counts.csv'
ace1.window.range <- paste(min(ace1.snps$Pos), max(ace1.snps$Pos), sep = '-')
system(paste('~/miniconda3/envs/gaard/bin/python extract_zarr_region.py', 
                                                 '~/scratch/VObs_GAARD',
                                                 study.id,
                                                 '2R',
                                                 ace1.window.range, 
                                                 ace1.allele.counts.filename,
                                                 sample.names.string, 
                                                 '../../data/tanzania.samples.meta.csv', 
                                                 'sample_id', 
                                                 'partner_sample_id', 
                                                 sep = ' '))
ace1.allele.counts <- fread(ace1.allele.counts.filename)
colnames(ace1.allele.counts)[1] <- 'pos'
ace1.allele.counts[, ID := paste('2R', pos, sep = ':')]
setkey(ace1.allele.counts, ID)

ace1.snps.ac <- ace1.allele.counts[ace1.snps$ID]
ace1.snps.ac[, chrom := '2R']
# as.numeric is needed because by taking the whole of the table's row as a vector, the numbers get turned into 
# characters
ace1.snps.ac[, alts := apply(.SD, 1, function(x) paste(x[..alt.allele.columns][as.numeric(x[..alt.counts.columns]) > 0], collapse = ','))]
ace1.snps.ac[, major.alt := apply(.SD, 1, function(x) paste(x[..alt.allele.columns][which.max(x[..alt.counts.columns])]))]
ace1.snps.ac[, info := ifelse(apply(.SD, 1, function(x) sum(x > 0)) > 1, 'Multiallelic', ''), .SDcols = alt.counts.columns]
ace1.vcf <- ace1.snps.ac[, .(chrom, pos, '.', ref = allele0, alt = alts, '.', '.', info)]
ace1.vcf.file <- 'temp.vcf'
ace1.vcf.output.file <- 'ace1.vcf'
write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO', file = vcf.file)
fwrite(ace1.vcf, ace1.vcf.file, sep = '\t', append = T, quote = F)

# Now run SNPEff on that. 
snp.eff.command <- paste(
	'java -jar ~/software/snpEff/snpEff.jar eff Anopheles_gambiae', 
	ace1.vcf.file, 
	'>', 
	ace1.vcf.output.file
)
cat('Running snpEff on Ace1 SNPs, using command', snp.eff.command, '.\n')
system(snp.eff.command)
system('rm temp.vcf')

ace1.snp.eff <- fread(ace1.vcf.output.file) %>%
                setnames('#CHROM', 'CHROM') %>%
				setnames(names(.), str_to_title(names(.)))

ace1.nonsyn <- merge(ace1.snps[, .(Chrom, Pos, MAF, P = logregP, Effect.direction = logreg.direction)], 
                     ace1.snp.eff[, .(Chrom, Pos, Info)], by = c('Chrom', 'Pos'))%>%
               .[grepl('missense', Info)] %>%
               .[, Effect.direction := as.factor(c('Resistant', 'Null', 'Susceptible')[Effect.direction + 2])] %>%
               .[, Nucleotide.change := sub('\\|.*', '', sub('.*\\|c\\.', '', Info))] %>%
               .[, Amino.acid.change := sub('\\|.*', '', sub('.*\\|p\\.', '', Info))] %>%
               .[, Info := NULL] %>%
               setnames(names(.), sub('\\.', ' ', names(.)))

fwrite(ace1.nonsyn, 'Ace1_nonsynonymous_SNPs.csv', sep = '\t')
