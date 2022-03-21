#-----------------------------------------------------------------------------
#
# Takes underlying files from QCEW on county-level income 
# and converts to panel data
#
#-----------------------------------------------------------------------------

library(stringr)
library(readxl)

# pick years from available data here
years <- 2000:2007


# create panel data
min_year <- min(years)

county_income_data <- data.frame()

for (year in years) {
  short_year <- str_sub(as.character(year), start=-2)
  file <- paste0("qcew_income_data/", year, "_all_county_high_level/allhlcn", short_year, ".xlsx")
  this_data <- read_excel(file)
  this_data <- this_data[this_data[,"Area Type"]=="County",] # drop national and state level dat
  this_data <- subset(this_data, Ownership=="Total Covered") # drop disaggregated industries
  this_data <- this_data[,c("St","Cnty","St Name","Area","Annual Average Pay")]
  this_data$year <- year
  colnames(this_data) <- c("state_code", "county_code", "state", "area", "annual_avg_pay", "year")
  if (year==min_year) {
    county_income_data <- this_data
  } else {
    county_income_data <- rbind.data.frame(county_income_data, this_data)
  }
}

county_income_data <- county_income_data[order(county_income_data$state, county_income_data$area, county_income_data$year),]

save(county_income_data, file="county_income_data.RData")
