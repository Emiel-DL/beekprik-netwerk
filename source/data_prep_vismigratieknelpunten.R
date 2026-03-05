library(sf)
library(tidyverse)
library(here)
library(mapview)

# inlezen knelpuntenlaag

knelpunten <- st_read(here("data", "ruw", "vismigratieknelpunten", "Export van knelpunten.shp"), quiet = T) %>%
  mutate(datum_vaststelling = dmy_hm(datum_vasts),
         datum_laatst_gewijzigd = dmy_hm(knp_laatst_)) %>%
  select(objectid, Knelpunt_nu, Status, datum_laatst_gewijzigd, knelpunttyp)

mapview(knelpunten)
