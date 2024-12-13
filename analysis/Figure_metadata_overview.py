import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import cm
meta = pd.read_csv('../data/ShanghaiDogsMetadata/SH_Dog_metadata_red.csv')

in2cm = 2.54
fig, ax = plt.subplots(figsize=(9/in2cm, 6/in2cm))
sex_dist = meta['Sex'].value_counts()

ax.barh(3, sex_dist['Female'], color=cm.Dark2.colors[0],)
ax.barh(3, sex_dist['Male'], color=cm.Dark2.colors[1], left=sex_dist['Female'])

size_dist = meta['Size_class'].value_counts()

ax.barh(2, size_dist['small'], color=cm.Dark2.colors[0],)
ax.barh(2, size_dist['medium'], color=cm.Dark2.colors[1], left=size_dist['small'])
ax.barh(2, size_dist['large'], color=cm.Dark2.colors[2], left=size_dist['small']+size_dist['medium'])

age_dist = meta['Age_classification'].replace(
        {
         'Mature_adult': 'Adult',
         'Early_senior': 'Senior',
         'Puppy': 'Juvenile',
         'Young_adult': 'Adult',
         'Geriatric': 'Senior',
         'Late_senior': 'Senior',
         }).value_counts()

ax.barh(1, age_dist['Juvenile'], color=cm.Dark2.colors[0],)
ax.barh(1, age_dist['Adult'], color=cm.Dark2.colors[1], left=age_dist['Juvenile'])
ax.barh(1, age_dist['Senior'], color=cm.Dark2.colors[2], left=age_dist['Juvenile']+age_dist['Adult'])
ax.barh(1, age_dist['unknown'], color=cm.Dark2.colors[7], left=age_dist['Juvenile']+age_dist['Adult']+age_dist['Senior'])

ax.set_ylim(-4, 4)
fig.tight_layout()
sns.despine()
fig.savefig('figures/dog_metadata.svg', dpi=300)

s = meta['Household.1'].value_counts().value_counts()
s = pd.DataFrame(
        {'Household size': s.index,
         'Number of households': s.values,
         'Number of dogs': s.values * s.index
         }
        ).set_index('Household size')

print(f'''Household size distribution:
{s.sort_index()}''')


