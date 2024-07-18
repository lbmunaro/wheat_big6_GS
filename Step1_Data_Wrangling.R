# Single trial analysis ----

# Objective ----
# - Organize and check phenotypic data

rm(list=objects()) # clean workspace

# Packages ----
library(tidyverse) # R packages for data science
library(janitor) # Simple Tools for Examining and Cleaning Dirty Data

# Functions ----
# Convert Yield to bu/acre
convYld<- function(y){
  x<- y/c(60 * 0.453592 * 2.47105)
  return(x)
}

# Convert test weight to lbs/bu
convTwt<- function(y){
  x<- y/1000 *2.2046 *35.2391
  return(x)
} 

# Phenotypic data ----

## Load and organize ----

### Wide format ----
YTpheno_raw_w <- #raw data wide format
  read.csv('Data/2024-07-12T214129phenotype_download.csv') |> # read csv file
  clean_names() |> # clean names
  dplyr::filter(study_name!="Big6_Scb_23") |> #remove scab nursery
  remove_empty(which = c('cols')) |> #remove columns entirely empty
  dplyr::select(!contains("growth_stage")) |> #remove columns that contains "growth_stage"
  mutate(trial=paste(study_year,location_name,sep = ' '),
         grain_yield_bu_ac=convYld(grain_yield_kg_ha_co_321_0001218), #convert grain yield to bu/ac
         test_weight_lb_bu=convTwt(grain_test_weight_g_l_co_321_0001210), #convert test weight to lb/bu
         heading_time_jd=heading_time_julian_date_jd_co_321_0001233, # shorter name
         plant_height_in=plant_height_cm_co_321_0001301/2.54, #convert plant_height to in
         maturity_jd = maturity_time_spike_estimation_julian_date_jd_co_321_0501101) |> # shorten name
  rename(year=study_year, location=location_name, study=study_name,germplasm=germplasm_name,
         observation=observation_unit_name, block=block_number, plot=plot_number,
         row=row_number, col=col_number) |> # shorten name
  dplyr::select(year, location, study, trial, germplasm, observation, replicate,
                block, plot, row, col, grain_yield_bu_ac:maturity_jd) |> #select columns
  arrange(study,row,col) |> #arrange by study, row, and col
  mutate_at(vars(year:col),as.factor) |> #convert to factor
  mutate_at(vars(grain_yield_bu_ac:maturity_jd),as.numeric) |> #convert to numeric
  # replace 0 by NA for response variables
  mutate_at(vars(grain_yield_bu_ac:maturity_jd),
            ~ifelse(.<=0,NA,.)) |>
  glimpse()

YTpheno_w <- YTpheno_raw_w |>
  # replace errors by NA
  mutate(grain_yield_bu_ac = ifelse(observation %in% c('Big6_Prn_23-257', 'Big6_Prn_23-32', 'Big6_Prn_23-333', 'Big6_Fre_23-60'), NA, grain_yield_bu_ac),
         test_weight_lb_bu = ifelse(observation %in% c('Big6_Prn_23-257', 'Big6_Prn_23-32', 'Big6_Prn_23-333', 'Big6_Prn_23-365', 'Big6_Prn_24-357', 'Big6_Prn_24-403'), NA, test_weight_lb_bu),
         plant_height_in = ifelse(observation == 'Big6_Urb_23-77', NA, plant_height_in)) |>
  # add row 56 & 57 col 34 from Big6_Urb_23.
  bind_rows(data.frame(year=as.factor(rep(2023,2)),
                       location=as.factor(rep('Urbana, IL',2)),
                       study=as.factor(rep('Big6_Urb_23',2)),
                       trial = as.factor(rep('2023 Urbana, IL',2)),
                       observation=as.factor(c('Big6_Urb_23-461','Big6_Urb_23-462')),
                       plot=as.factor(c(461,462)),
                       row=as.factor(c(56,57)),
                       col=as.factor(c(34,34)))) |>
  arrange(study,row,col) |> #arrange by study, row, and col
  glimpse()

### Long format ----
YTpheno_l <- YTpheno_w |>
  # Transform from wide to long format
  pivot_longer(cols = grain_yield_bu_ac:maturity_jd,
               names_to = 'trait', values_to = 'value', values_drop_na = F) |>
  mutate(trait=as.factor(trait)) |> # trait as factor
  arrange(study, trait, as.numeric(row), as.numeric(col)) |> # arrange
  group_by(study, trait) |>
  filter(!is.na(mean(value, na.rm = TRUE))) |> # filter out traits without pheno obs.
  ungroup() |>
  glimpse()

### Data visualization ----

#### Boxplot ----
YTpheno_l |>
  ggplot(aes(x=study, y=value)) +
  geom_boxplot() +
  facet_wrap(~trait, scales = 'free', ncol = 1)

#### Histogram ----
YTpheno_l |>
  ggplot(aes(x=value)) +
  geom_histogram() +
  facet_wrap(~trait:study, scales = 'free')

#### Tile plot ----

#### Block number ----
YTpheno_l |>
  filter(trait=='grain_yield_bu_ac') |>
  ggplot(aes(x=col,y=row, fill=block)) +
  geom_tile() +
  facet_wrap(~study, ncol = 2) +
  theme_bw()

#### Grain yield ----
YTpheno_l |>
  filter(trait=='grain_yield_bu_ac') |>
  ggplot(aes(x=col,y=row, fill=value)) +
  geom_tile() +
  facet_wrap(~study, ncol = 2, scales = 'free') +
  scale_fill_viridis_c(name= 'Grain yield',
                       option = 'viridis',
                       na.value='white') +
  theme_bw()

#### Test weight ----
YTpheno_l |>
  filter(trait=='test_weight_lb_bu') |>
  ggplot(aes(x=col,y=row, fill=value)) +
  geom_tile() +
  facet_wrap(~study, ncol = 2, scales = 'free') +
  scale_fill_viridis_c(name= 'Test weight',
                       option = 'viridis',
                       na.value='white') +
  theme_bw()

#### Heading time ----
YTpheno_l |>
  filter(trait=='heading_time_jd') |>
  ggplot(aes(x=col,y=row, fill=value)) +
  geom_tile() +
  facet_wrap(~study, ncol = 2, scales = 'free') +
  scale_fill_viridis_c(name= 'Heading time',
                       option = 'viridis',
                       na.value='white') +
  theme_bw()

#### Plant height ----
YTpheno_l |>
  filter(trait=='plant_height_in') |>
  ggplot(aes(x=col,y=row, fill=value)) +
  geom_tile() +
  facet_wrap(~study, ncol = 2, scales = 'free') +
  scale_fill_viridis_c(name= 'Plant height',
                       option = 'viridis',
                       na.value='white') +
  theme_bw()

#### Maturity ----
YTpheno_l |>
  filter(trait=='maturity_jd') |>
  ggplot(aes(x=col,y=row, fill=value)) +
  geom_tile() +
  facet_wrap(~study, ncol = 2, scales = 'free') +
  scale_fill_viridis_c(name= 'Maturity',
                       option = 'viridis',
                       na.value='white') +
  theme_bw()


# Map of the locations ----
unique(YTpheno_w$location)

library(maps)
library(mapdata)
library(ggrepel)

locations <- data.frame(
  city = c("Frankenmuth", "Fremont", "Ithaca", "West Lafayette", "Mason",
           "Neoga", "Princeton", "Urbana", "Vincennes", "Wooster"),
  state = c("MI", "OH", "NY", "IN", "MI", "IL", "KY", "IL", "IN", "OH"),
  lat = c(43.3317, 41.3506, 42.4430, 40.4259, 42.5795,
          39.4437, 37.8644, 40.1106, 38.6773, 40.8051),
  lon = c(-83.7372, -83.1138, -76.5019, -86.9081, -84.4433,
          -88.4528, -87.8669, -88.2073, -87.5286, -81.9381))

# Get the map data for the US
usa <- map_data("state")

# Filter the map data to include only the specified states
states_to_include <- c("new york", "ohio", "michigan", "indiana", "illinois", "kentucky")
filtered_map_data <- usa %>% filter(region %in% states_to_include)

# Create the map plot
ggplot() +
  geom_polygon(data = filtered_map_data, aes(x = long, y = lat, group = group), fill = "lightgray", color = "white") +
  geom_point(data = locations, aes(x = lon, y = lat), color = "#FF5F0F", size = 1) +
  geom_text_repel(data = locations, aes(x = lon, y = lat, label = city), size = 1.5, nudge_y = 0.2) +
  coord_fixed(1.3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x = "Longitude",y = "Latitude")
ggsave('Figures/locations_map.png') # save plot

# Save image
save.image('Data/Big6_pheno_data.RData')