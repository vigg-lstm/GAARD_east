from sys import argv
import numpy as np
import pandas as pd
import malariagen_data
import allel

if (len(argv) == 4):
    cohort = argv[1]
    contig = argv[2]
    window_size = argv[3]
else:
    raise Exception("Fail. There should be three command line arguments (cohort, contig, window_size).")

window_size = int(window_size)

ag3 = malariagen_data.Ag3(pre=True)
metadata = ag3.sample_metadata("3.7")
metadata = metadata.query("sex_call == 'F'")
sibs = pd.read_csv("../NGSrelate/full_relatedness_tanzania/sib_group_table.csv", sep="\t")
# Remove the sibs and the samples with odd relatedness
exclude = pd.concat([
    sibs.query("keep == False")['sample.name'],
    pd.read_csv('../NGSrelate/full_relatedness_tanzania/samples_to_exclude.csv', header = None)[0]
])

metadata = metadata.query("partner_sample_id not in @exclude").reset_index(drop=True)

# Identify the samples from the correct cohort
phenotypes = pd.read_csv('../data/sample_phenotypes_EA.csv', sep = '\t', na_filter = False)
phenotypes['population'] = phenotypes[['location', 'species', 'insecticide']].agg('.'.join, axis = 1)
phenotypes = phenotypes.loc[phenotypes['population'].str.contains(cohort), :]
metadata = metadata.query(f"partner_sample_id in {phenotypes['specimen'].to_list()}")
cohort_ids = metadata['partner_sample_id'].to_list()

pos, h12 = ag3.h12_gwss(
    contig=contig,
    analysis="arab",
    window_size=window_size,
    sample_sets="1246-VO-TZ-KABULA-VMF00185",
    sample_query=f"partner_sample_id in {cohort_ids}"
)

h12_df = pd.DataFrame({'midpoint': pos,
                       'h12': np.round(h12, 3)})

h12_df.to_csv(f"H12_{cohort}_{contig}.csv", sep="\t", index = False)
