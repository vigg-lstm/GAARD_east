library(data.table)
library(magrittr)

combine.haplotypes <- function(sample.set, hap.files){
	haplotype.vcfs <- grep(paste(sample.set, '.*\\.vcf', sep = ''), hap.files, value = T) 
	haplotype.csvs <- grep(paste(sample.set, '.*hap_sequence\\.csv', sep = ''), hap.files, value = T)
	allsites.vcf <- lapply(haplotype.vcfs, fread, select = c('#CHROM', 'POS', 'INFO')) %>%
	                rbindlist() %>%
	                setnames(c('#CHROM', 'POS', 'INFO'), c('Chrom', 'Pos', 'Info'))
	allsites.csv <- lapply(haplotype.csvs, fread) %>%
	                rbindlist(fill = T)
	full.table <- merge(allsites.csv, allsites.vcf)
	fwrite(full.table, file = paste(sample.set, '_haplotype_SNPs.csv', sep = ''), sep = '\t')
}

sample.sets <- c('Moshi_arabiensis_Delta', 
                 'Moshi_arabiensis_PM', 
                 'Muleba_arabiensis_Delta')

sapply(sample.sets, 
       combine.haplotypes, 
       hap.files = list.files('../../haplotypes/', full.names = T)
)


