print("Loading packages...")
# List of required packages
packages <- c("stringdist", 
              "tidyverse", 
              "AER", 
              "sandwich", 
              "lmtest", 
              "plm", 
              "readr", 
              "geosphere")
# Loop through the list and install any package that isn't installed yet, then load it
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}
print("packages loaded successfully")

print("Loading data...")
# Read NDR, Geocoded and DOE data
SA <- read_csv("Data/SA_Data_all.csv") # Replace by NDR
Geocoded_Station_State <- read_csv("Data/Geocoded_Station_State.csv")
dt <- read_csv("Data/data_wind_com.csv")

print("Data loaded successfully!")

### 01. Data Cleaning

# Cleaning for NDR
relevant_columns <- c("PATIENTID", "R_date", "D_date", "status", "Region",
                      "State", "District",  "Age_At_DC", "Latitude", "Longitude")
SA_subset <- SA %>%
  select(all_of(relevant_columns))

SA_subset <- SA_subset %>% mutate(State = tolower(State))

SA_subset <- SA_subset %>%
  mutate(D_date = as.Date(D_date, format = "%d/%m/%Y"))

SA_subset <- SA_subset %>%
  mutate(Latitude = as.numeric(Latitude), Longitude = as.numeric(Longitude))

# Cleaning for Geo
Geocoded_Station_State <- Geocoded_Station_State %>%
  mutate(lat = as.numeric(lat), lon = as.numeric(lon))

# Cleaning for DOE
dt <- dt %>%
  mutate(date = as.Date(date))

# Extract month and year from the date
dt <- dt %>%
  mutate(
    year = year(date),
    month = month(date)
  )

# Convert month and year to factor with consistent levels
dt$year <- factor(dt$year, levels = unique(dt$year))
dt$month <- factor(dt$month, levels = unique(dt$month))

dt <- dt %>%
  filter(date >= as.Date("2017-07-04"))

# Check NAs
dt %>%
  summarise(across(everything(), 
                   ~sum(is.na(.)))) %>%
  pivot_longer(everything(), 
               names_to = "Variable", 
               values_to = "NA_Count") %>%
  mutate(NA_Percentage = (NA_Count/nrow(dt))*100) %>%
  arrange(desc(NA_Count))

dt %>%
  select(PM10AVG, SO2AVG, NO2AVG, NOAVG, NOXAVG, O3AVG, COAVG, PM2.5AVG) %>%
  summarise(across(everything(), 
                   ~sum(is.na(.)))) %>%
  pivot_longer(everything(), 
               names_to = "Pollutant", 
               values_to = "NA_Count") %>%
  mutate(NA_Percentage = (NA_Count/nrow(dt))*100) %>%
  arrange(desc(NA_Percentage))

### 02. Data Merging

print("Running the loop to match geodata...")
# Create a function to find the nearest station for each patient
find_nearest_station <- function(patient_lat, 
                                 patient_lon, 
                                 stations_lats, 
                                 stations_lons) {
  # Calculate distances between patient and all stations
  distances <- distGeo(cbind(stations_lons, 
                             stations_lats), 
                       c(patient_lon, 
                         patient_lat))
  # Return the index of the nearest station
  nearest_station_index <- which.min(distances)
  return(nearest_station_index)
}

# Apply the proximity matching to assign each patient to the nearest station
SA_subset <- SA_subset %>%
  rowwise() %>%
  mutate(
    nearest_station_index = find_nearest_station(Latitude, 
                                                 Longitude, 
                                                 Geocoded_Station_State$lat, 
                                                 Geocoded_Station_State$lon),
    nearest_station = Geocoded_Station_State$Location[nearest_station_index],  # Assign nearest station name
    nearest_station_lat = Geocoded_Station_State$lat[nearest_station_index],  # Assign nearest station latitude
    nearest_station_lon = Geocoded_Station_State$lon[nearest_station_index]   # Assign nearest station longitude
  )

print("Geodata matched successfully!")


# Convert nearest_station_index into a numeric type 
SA_subset <- SA_subset %>%
  mutate(nearest_station_index = as.numeric(nearest_station_index))

clean_location_names <- function(location) {
  location %>%
    tolower() %>%
    # Remove school/institution prefixes
    str_remove_all("(?i)sek\\.? ?(men|keb)\\.? ?") %>%
    str_remove_all("(?i)institut\\w* ?") %>%
    str_remove_all("(?i)kolej\\w* ?") %>%
    str_remove_all("(?i)universiti\\w* ?") %>%
    # Remove common location prefixes
    str_remove_all("(?i)pejabat\\w* ?") %>%
    str_remove_all("(?i)kuarters\\w* ?") %>%
    str_remove_all("(?i)komplek\\w* ?") %>%
    # Keep main location name
    str_extract("([\\w\\s]+),\\s*[\\w\\.\\s]+$") %>%
    trimws()
}

# Get unique locations first
unique_dt_locations <- unique(dt$LOCATION)
unique_sa_locations <- unique(SA_subset$nearest_station)

# Modified match_locations function
match_locations <- function(dt_locations, sa_locations) {
  # Clean both location sets
  dt_clean <- sapply(dt_locations, clean_location_names)
  sa_clean <- sapply(sa_locations, clean_location_names)
  
  # Create distance matrix
  distances <- stringdistmatrix(dt_clean, sa_clean, method = "jw")
  
  # Find best matches
  matches <- data.frame(
    dt_location = dt_locations,
    dt_clean = dt_clean,
    sa_location = sa_locations[apply(distances, 1, which.min)],
    sa_clean = sa_clean[apply(distances, 1, which.min)],
    distance = apply(distances, 1, min)
  )
  
  return(matches %>% arrange(distance))
}

# Run with unique locations
matches <- match_locations(unique_dt_locations, unique_sa_locations)

# Create mapping dictionary
location_mapping <- c(
  "seberang jaya, pulau pinang" = "seberang perai, penang",
  "balik pulau, pulau pinang" = "balik pulau, penang",
  "minden, pulau pinang" = "minden, penang"
  # Add other mappings as needed
)

# Apply mapping and merge
dt <- dt %>%
  mutate(LOCATION_clean = tolower(LOCATION))

SA_subset <- SA_subset %>%
  mutate(nearest_station_clean = tolower(nearest_station))

# Update columns using mapping
for(old_name in names(location_mapping)) {
  dt$LOCATION_clean[dt$LOCATION_clean == old_name] <- location_mapping[old_name]
  SA_subset$nearest_station_clean[SA_subset$nearest_station_clean == old_name] <- location_mapping[old_name]
}

# Find common stations
common_stations <- intersect(unique(dt$LOCATION_clean), 
                             unique(SA_subset$nearest_station_clean))

# Filter both datasets for common stations
SA_subset_common <- SA_subset %>%
  filter(nearest_station_clean %in% common_stations)

dt_common <- dt %>%
  filter(LOCATION_clean %in% common_stations)

# Merge datasets
merged_death_data <- SA_subset_common %>%
  left_join(dt_common, 
            by = c("nearest_station_clean" = "LOCATION_clean",
                   "D_date" = "date"))

merged_death_data <- merged_death_data %>%
  filter(D_date >= as.Date("2017-07-04"))

# Create daily death counts from D_date
daily_deaths <- SA_subset %>%
  count(D_date, nearest_station_clean) %>%
  rename(date = D_date, deaths = n)

# Merge with environmental data
merged_daily <- dt_common %>%
  left_join(daily_deaths, 
            by = c("date", 
                   "LOCATION_clean" = "nearest_station_clean")) %>%
  mutate(deaths = replace_na(deaths, 0))  # Days without deaths get 0

print("Merge performed successfully!")

write.csv(merged_daily, 'Data/Mid_process_data/merged_death_data_from_script_01.csv')
