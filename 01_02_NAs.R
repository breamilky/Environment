# Function to read, clean, and ensure consistent data types
read_and_clean <- function(sheet) {
  df <- read_excel(single_path, sheet = sheet)
  # Convert 'NA' strings to actual NA, handle date columns separately
  df <- df %>%
    mutate(across(everything(), as.character)) %>%
    mutate(across(everything(), ~na_if(.x, "NA")))
  return(df)
}
# Read all sheets and combine them into one dataframe
NA_combined_df <- bind_rows(lapply(single_sheet, read_and_clean))

# Convert necessary columns back to their appropriate types
numeric_columns <- c("PM10 AVG (µg/m3)", "PM2.5 AVG (µg/m3)", "SO2 AVG (ppm)", 
                     "NO2 AVG (ppm)", "NO AVG (ppm)", "NOX AVG (ppm)", 
                     "O3 AVG (ppm)", "CO AVG (ppm)", "WIND DIRECTION AVG(°)", "WIND SPEED AVG(m/s)")

datetime_columns <- c("DATE TIME")

NA_combined_df <- NA_combined_df %>%
  mutate(across(all_of(numeric_columns), ~ as.numeric(.x))) %>%
  mutate(across(all_of(datetime_columns), ~ as.POSIXct(.x, format = "%d/%m/%Y")))
na_counts <- sapply(NA_combined_df, function(x) sum(is.na(x)))

# Print the NA counts
print(na_counts)