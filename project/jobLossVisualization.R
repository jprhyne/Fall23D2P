library(dplyr)
library(sf)
library(readr)
# County population data manually entered from 2022
load("pop.rda")
# Read in the provided geojson file
USMap = st_read("all_jobs_sum_job_loss_county.geojson")
# Filter out to only the Colorado level
COMap = USMap %>% filter(state_name=="Colorado")
# Frees the memory used by USMap above. Helpful for lower RAM devices
USMap = ""
# Renames the nondescriptive variable to something that we can tell what it means
names(COMap)[names(COMap) == "X000"] = "Total"
# Creates a new variable that will be the rates. For every region we divide by the population to get an estimate of the jobs lost per
# person
COMap$Rates=COMap$Total/pop
# Plot the exploratory data
plot(COMap["Rates"])
# Give a descriptive title
title(main="Rates of lost jobs")

save(COMap, file="COMap.rda")
