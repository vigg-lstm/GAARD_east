library(data.table)
library(magrittr)
library(stringr)

sequenced.meta <- fread('../../data/tanzania.samples.meta.csv')[sex_call == 'F'] %>%
                  setnames('partner_sample_id', 'sample.id')
phen <- fread('../../data/sample_phenotypes_EA.csv') %>%
              setnames(c('specimen', 'exposure_time'), c('sample.id', 'exposure.time'))
sibs <- fread('../../NGSrelate/full_relatedness_tanzania/sib_group_table.csv') %>%
              setnames('sample.name', 'sample.id')
accessions <- fread('../../data/wgs_lanelets_ena_tanzania.csv')

sequenced.meta[, collection.date := paste(str_pad(month, 2, pad = '0'), year, sep = '/')]
output.table <- merge(sequenced.meta[, .(sample.id, sample_id, country, location, latitude, longitude, collection.date)],
                      phen[, !c('location', 'country', 'plate')],
					  all = F, by = 'sample.id') %>%
				merge(., sibs[, .(sample.id, sib.group.id, exclude.as.sib = !keep)],
					  all = T, by = 'sample.id') %>%
				merge(., accessions[, .(sample_id, run.accession = run_ena)],
					  all = F, by = 'sample_id') %>%
				.[, !c('sample_id')]

output.table[is.na(exclude.as.sib), exclude.as.sib := F]
			
fwrite(output.table, 'metadata.csv', sep = '\t')

