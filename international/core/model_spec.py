import numpy as np

def z(a):
    return (a - a.mean(0)) / a.std(0, ddof=1)

def rcs(age):
    t = np.array([18.0, 30.0, 45.0, 60.0])
    x = np.asarray(age, float)

    def d(k):
        return np.maximum(x - k, 0) ** 3
    return np.column_stack([x - 35] + [(d(t[j]) - d(t[-2]) * (t[-1] - t[j]) / (t[-1] - t[-2]) + d(t[-1]) * (t[-2] - t[j]) / (t[-1] - t[-2])) / (t[-1] - t[0]) ** 2 for j in range(len(t) - 2)])

def individual_columns(d, nations):
    g = d['g']
    count = len(d['C'])
    columns = {'intercept': np.ones(len(g))}

    def pair(name, value):
        means = np.array([value[g == j].mean() for j in range(count)])[g]
        columns[name + '_within'] = z(value - means)
        columns[name + '_between'] = z(means)
    gmt = d['harm'].mean(1)
    for (name, value) in [('female', d['female']), ('GMT', gmt), ('TIV', d['tiv'])]:
        pair(name, value)
    basis = rcs(d['age'])
    for j in range(3):
        name = ('age_rcs' if nations else 'age') + str(j + 1)
        pair(name, basis[:, j])
        pair(name + ('_x_female' if nations else '_sex'), basis[:, j] * (d['female'] - 0.5))
    pair('age_GMT', z(d['age']) * z(gmt))
    return columns

NATIONS_DEMOGRAPHICS = ['intercept', 'female_within', 'GMT_within', 'TIV_within', 'age_rcs1_within', 'age_rcs2_within', 'age_rcs3_within', 'age_rcs1_x_female_within', 'age_rcs2_x_female_within', 'age_rcs3_x_female_within', 'age_GMT_within', 'female_between', 'TIV_between', 'age_rcs1_between', 'age_GMT_between', 'age_rcs2_between', 'age_rcs3_between']

USA_DEMOGRAPHICS = ['intercept', 'female_within', 'female_between', 'GMT_within', 'TIV_within', 'age1_within', 'age2_within', 'age3_within', 'age1_sex_within', 'age1_sex_between', 'age2_sex_within', 'age3_sex_within', 'age_GMT_within']

NATIONS_BLOCKS = {'socioeconomic': ['HDI2020', 'HDI1990', 'PPPrw'], 'climate': ['tmp', 'frs', 'wet'], 'land_resources': ['freshwater', 'arable_land', 'food_production'], 'geography': ['abslatitude'], 'PM2_5': ['PM2_5'], 'forest_cover': ['forest_cover']}

USA_MACROS = ['humidity', 'disaster_rate_per_people', 'land_forested', 'population_density', 'precipitation', 'air_quality_index', 'PPP_2017', 'highest_elevation_ft']
