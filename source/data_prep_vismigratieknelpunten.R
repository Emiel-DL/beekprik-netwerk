library(sf)
library(tidyverse)
library(here)
library(mapview)

# inlezen knelpuntenlaag

knelpunten <- st_read(here("data", "ruw", "vismigratieknelpunten", "Export van knelpunten.shp"), quiet = T) %>%
  janitor::clean_names() %>%
  mutate(datum_vaststelling = dmy_hm(datum_vasts),
         datum_laatst_gewijzigd = dmy_hm(knp_laatst)) %>%
  select(objectid, knelpunt_nu, status, datum_laatst_gewijzigd, knelpunttyp, vhag, vhas, waterlichaa) %>%
  rename(knelpunt_nummer = knelpunt_nu,
         type = knelpunttyp,
         owl = waterlichaa) %>%
  st_transform(., crs = 31370)

library(dplyr)
library(tidyr)
library(lubridate)

# 1. Definieer het bereik
startjaar_analyse <- min(year(knelpunten$datum_laatst_gewijzigd), na.rm = TRUE) # Pas dit aan naar wens
huidig_jaar <- year(today())

knelpunten_tijdreeks <- knelpunten %>%
  # Haal het jaar uit de wijzigingsdatum
  mutate(jaar_laatst_gewijzigd = lubridate::year(datum_laatst_gewijzigd)) %>%
  # Bepaal tot welk jaar het knelpunt bestond/bestaat
  mutate(eindjaar = case_when(
    status == "Opgelost" ~ jaar_laatst_gewijzigd - 1,
    TRUE ~ huidig_jaar
  )) %>%
  # Creëer een lijst van jaren voor elk object
  rowwise() %>%
  mutate(jaar = list(seq(startjaar_analyse, eindjaar))) %>%
  ungroup() %>%
  # Vouw de lijst uit naar aparte rijen
  unnest(jaar) %>%
  # Optioneel: sorteer op jaar en id
  arrange(jaar, knelpunt_nummer)

# Check het resultaat voor één specifiek knelpunt
knelpunten_tijdreeks %>% filter(knelpunt_nummer == "5003-010") %>% View

head(knelpunten_tijdreeks)
