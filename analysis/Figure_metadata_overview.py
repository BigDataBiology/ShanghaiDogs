import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm

plt.rcParams['svg.fonttype'] = 'none' #to avoid transforming the font to plot

meta = pd.read_csv('data/ShanghaiDogsMetadata/SH_Dog_metadata_red.csv')

in2cm = 2.54
fig, ax = plt.subplots(figsize=(7.5/in2cm, 6/in2cm))

# Sex
sex_dist = meta['Sex'].value_counts()
ax.barh(3, sex_dist['Female'], color=cm.Dark2.colors[0],)
ax.barh(3, sex_dist['Male'], color=cm.Dark2.colors[1], left=sex_dist['Female'])

ax.text(sex_dist['Female'] / 2, 3, "female", va='center', ha='center', color='white', fontsize=9)
ax.text(sex_dist['Female'] + sex_dist['Male'] / 2, 3, "male", va='center', ha='center', color='white', fontsize=9)

# Size
size_dist = meta['Size_class'].value_counts()
ax.barh(2, size_dist['small'], color=cm.Dark2.colors[0],)
ax.barh(2, size_dist['medium'], color=cm.Dark2.colors[1], left=size_dist['small'])
ax.barh(2, size_dist['large'], color=cm.Dark2.colors[2], left=size_dist['small']+size_dist['medium'])

ax.text(size_dist['small'] / 2, 2, "S", va='center', ha='center', color='white', fontsize=9)
ax.text(size_dist['small'] + size_dist['medium'] / 2, 2, "M", va='center', ha='center', color='white', fontsize=9)
ax.text(size_dist['small'] + size_dist['medium'] + size_dist['large'] / 2, 2, "L", va='center', ha='center', color='white', fontsize=9)

# Age
age_dist = meta['Animal_age_simplified'].value_counts()
ax.barh(1, age_dist['Young'], color=cm.Dark2.colors[0],)
ax.barh(1, age_dist['Adult'], color=cm.Dark2.colors[1], left=age_dist['Young'])
ax.barh(1, age_dist['Senior'], color=cm.Dark2.colors[2], left=age_dist['Young']+age_dist['Adult'])
ax.barh(1, age_dist['unknown'], color=cm.Dark2.colors[7], left=age_dist['Young']+age_dist['Adult']+age_dist['Senior'])

ax.text(age_dist['Young'] / 2, 1, "young", va='center', ha='center', color='white', fontsize=9)
ax.text(age_dist['Young'] + age_dist['Adult'] / 2, 1, "adult", va='center', ha='center', color='white', fontsize=9)
ax.text(age_dist['Young'] + age_dist['Adult'] + age_dist['Senior'] / 2, 1, "senior", va='center', ha='center', color='white', fontsize=9)

# Format
ax.set_ylim(-4, 4)
ax.axis('off')
fig.tight_layout()
#plt.show()
fig.savefig('intermediate-outputs/figures/dog_metadata.svg', dpi=300)

s = meta['Household.1'].value_counts().value_counts()
s = pd.DataFrame(
        {'Household size': s.index,
         'Number of households': s.values,
         'Number of dogs': s.values * s.index
         }
        ).set_index('Household size')

print(f'''Household size distribution:
{s.sort_index()}''')


