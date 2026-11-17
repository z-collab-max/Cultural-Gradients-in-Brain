"""Raw macro fields, used only to assemble PCA inputs."""
CHINA_CITY = ['absolute_latitude', 'terrain_ruggedness', 'historical_humidity',
              'historical_temperature', 'historical_sunshine', 'historical_precipitation',
              'population_density', 'green_space', 'natural_disaster_frequency_per_capita']
CHINA_PROVINCE = ['hdi_education', 'hdi_income']
NATIONS_BLOCKS = {'socioeconomic': ['HDI2020', 'HDI1990', 'PPPrw'],
                  'climate': ['tmp', 'frs', 'wet'],
                  'land_resources': ['freshwater', 'arable_land', 'food_production'],
                  'geography': ['abslatitude'], 'PM2_5': ['PM2_5'], 'forest_cover': ['forest_cover']}
USA_MACROS = ['humidity', 'disaster_rate_per_people', 'land_forested', 'population_density',
              'precipitation', 'air_quality_index', 'PPP_2017', 'highest_elevation_ft']
