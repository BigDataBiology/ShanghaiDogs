import numpy as np
import pandas as pd
from glob import glob
from matplotlib import pyplot as plt
import seaborn as sns

shd_mapstats = glob('outputs/*_shd_mapstats.txt')

fraction = {}
for file in sorted(shd_mapstats):
    sample = file.split('_')[0].split('/')[-1]
    fraction[sample] = pd.read_csv(file, sep='\t', index_col=0).T.eval('aligned/total').item()


fraction = pd.Series(fraction)
fraction = fraction.sort_values()
pd.DataFrame({'fraction': fraction}) \
    .reset_index() \
    .rename(columns={'index': 'sample'}) \
    .to_csv('outputs/fraction_shd_mapped.tsv', sep='\t', index=False)

