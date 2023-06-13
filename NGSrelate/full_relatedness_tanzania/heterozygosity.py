import pandas as pd
import numpy as np
import malariagen_data
import allel
import xarray as xr
#import seaborn as sns
#import matplotlib as plt
from scipy.stats import pearsonr

ag3 = malariagen_data.Ag3(pre = True)
ag3

meta = ag3.sample_metadata(sample_sets = '1246-VO-TZ-KABULA-VMF00185')
meta.set_index('sample_id', inplace = True)

snp_calls = (
    ag3.snp_calls(region = '3L',
                  sample_sets = '1246-VO-TZ-KABULA-VMF00185')
)

sample_ids = snp_calls.sample_id.values
sample_names = meta.loc[sample_ids, 'partner_sample_id'].values

snp_genotypes = allel.GenotypeArray(snp_calls.call_genotype.values)
snp_genotypes

missing_filter = snp_genotypes.is_missing().any(1)
snp_genotypes_filtered = snp_genotypes[~missing_filter, :, :]
snp_genotypes_filtered
hets = snp_genotypes_filtered.is_het()
heterozygosities = np.mean(hets, 0)

hets_df = pd.DataFrame({'het': heterozygosities}, index = sample_names)
hets_df.to_csv('heterozygosities.csv')
