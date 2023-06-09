from sys import argv
import numpy as np
import pandas as pd
import malariagen_data
import allel

if (len(argv) == 5):
    cohort1 = argv[1]
    cohort2 = argv[2]
    contig = argv[3]
    window_size = argv[4]
else:
    raise Exception("Fail. There should be four command line arguments (cohort1, cohort2, contig, window_size).")

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
phen1 = phenotypes.query(f"population == @cohort1")
meta1 = metadata.query(f"partner_sample_id in {phen1['specimen'].to_list()}")
cohort1_ids = meta1['partner_sample_id'].to_list()
phen2 = phenotypes.query(f"population == @cohort2")
meta2 = metadata.query(f"partner_sample_id in {phen2['specimen'].to_list()}")
cohort2_ids = meta2['partner_sample_id'].to_list()

pos, h1x = ag3.h1x_gwss(
    contig=contig,
    analysis="arab",
    window_size=window_size,
    sample_sets="1246-VO-TZ-KABULA-VMF00185",
    cohort1_query=f"partner_sample_id in {cohort1_ids}",
    cohort2_query=f"partner_sample_id in {cohort2_ids}"
)

h1x_df = pd.DataFrame({'midpoint': pos,
                       'h1x': np.round(h1x, 3)})

h1x_df.to_csv(f"H1x_{cohort1}_vs_{cohort2}_{contig}.csv", sep="\t", index = False)
