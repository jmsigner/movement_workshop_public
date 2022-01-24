library(tidyverse)
library(amt)

# Elephant data ------

URL <- "https://www.datarepository.movebank.org/bitstream/handle/10255/move.373/Elliptical%20Time-Density%20Model%20%28Wall%20et%20al.%202014%29%20African%20Elephant%20Dataset%20%28Source-Save%20the%20Elephants%29.csv"

elephants <- read_csv(url(URL))
elephants

# Time zone GMT
# UTM zone: 30N
elephants %>% select(long = `location-long`, lat = `location-lat`, 
                     temperature = `external-temperature`, 
                     timestamp = `study-local-timestamp`,
                     id = `individual-local-identifier`) %>% 
  write_csv("data/elephants.csv")

# Cougar and coyote data ----
Mahoney_URL <- "https://www.datarepository.movebank.org/bitstream/handle/10255/move.590/Site%20fidelity%20in%20cougars%20and%20coyotes%2c%20Utah_Idaho%20USA%20%28data%20from%20Mahoney%20et%20al.%202016%29-gps.csv"

coy_coug <- read_csv(url(Mahoney_URL))

# Time zone US/Mountain
# UTM zone: 12N
s_ut <- coy_coug %>% 
  select(long = `location-long`,
         lat = `location-lat`,
         t = timestamp,
         species = `individual-taxon-canonical-name`,
         id = `individual-local-identifier`) %>%
  # Subset just to southern Utah
  filter(lat < 41) %>%
  # Format as track
  make_track(long, lat, t, id = id, crs = 4326) %>% 
  # Convert to UTM
  transform_coords(32612)

write_csv(s_ut, "data/coyote_cougar.csv")

