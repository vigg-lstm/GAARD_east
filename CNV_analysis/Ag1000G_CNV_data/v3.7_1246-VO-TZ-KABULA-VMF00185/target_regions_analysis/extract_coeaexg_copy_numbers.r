library(abind)
library(magrittr)

load('target_regions_analysis.Rdata')

cnv.pos <- c(37282000, 37295000)

cnv.windows <- compact.hmm$coeaexg[[1]]$Position %>%
               {. >= cnv.pos[1] & . <= cnv.pos[2]}
coeaexg.cnv.cov <- lapply(compact.hmm$coeaexg, function(x) x[cnv.windows, c('Normalised_coverage', 'CNV')]) %>%
                   abind(rev.along = 0) %>%
                   aperm(c(1, 3, 2))
attributes(coeaexg.cnv.cov)$Position <- as.numeric(dimnames(coeaexg.cnv.cov)[[1]])

saveRDS(coeaexg.cnv.cov, file = 'Coeaexg_cnv_coverage.RDS')


