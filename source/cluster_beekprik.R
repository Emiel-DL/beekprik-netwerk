library(sf)
library(igraph)
library(tidyverse)
library(progress)
library(here)
library(mapview)
library(lwgeom)
library(mapview)

beekprik <- read.csv(here("data", "ruw", "beekprik", "Finale_data_beekprik_aan_en_afwezigheid_1990_2023.csv")) %>%
  janitor::clean_names() %>%
  st_as_sf(., coords = c("lambert_x", "lambert_y"), crs = 31370)

beekprik_aanwezig <- beekprik %>%
  filter(aanwezig == 1)

mapview(beekprik_aanwezig)

# netwerk bouwen

fd <- st_read(here("data", "ruw", "netwerk", "Flow_direction_coordinates.shp"), quiet = T)
nodes <- st_read(here("data", "ruw", "netwerk", "vha_network_junctions.shp"), quiet = T)

split_lines_equal_length <- function(lines_sf, max_length = 200) {

  message(paste("Lijnen splitsen in segmenten van max", max_length, "meter..."))

  # 1. Voorbereiding
  lines_sf <- st_cast(lines_sf, "LINESTRING", warn = FALSE)
  lines_sf$orig_len <- as.numeric(st_length(lines_sf))
  lines_sf$orig_id  <- 1:nrow(lines_sf)

  # 2. Scheid lijnen
  short_lines <- lines_sf %>% filter(orig_len <= max_length)
  long_lines  <- lines_sf %>% filter(orig_len > max_length)

  if (nrow(long_lines) == 0) {
    return(lines_sf)
  }

  message(paste(" ->", nrow(long_lines), "lange lijnen worden opgeknipt..."))

  # 3. Bereken aantal segmenten en blaas tabel op
  long_lines$n_segments <- ceiling(long_lines$orig_len / max_length)

  expanded_lines <- long_lines[rep(row.names(long_lines), long_lines$n_segments), ]

  expanded_lines <- expanded_lines %>%
    group_by(orig_id) %>%
    mutate(seg_num = row_number()) %>%
    ungroup()

  # 4. Bereken fracties
  expanded_lines <- expanded_lines %>%
    mutate(
      step_size = max_length / orig_len,
      start_frac = (seg_num - 1) * step_size,
      end_frac   = seg_num * step_size
    )

  # Correctie voor afrondingsfouten (nooit meer dan 100% van de lijn)
  expanded_lines$end_frac[expanded_lines$end_frac > 1] <- 1

  message(" -> Geometrieën daadwerkelijk knippen (dit kan even duren)...")

  list_geoms <- mapply(
    function(geom, f_start, f_end) {
      st_linesubstring(geom, from = f_start, to = f_end, tolerance = 0.001)
    },
    expanded_lines$geometry,   # Input 1: De geometrie
    expanded_lines$start_frac, # Input 2: Start
    expanded_lines$end_frac,   # Input 3: Eind
    SIMPLIFY = FALSE           # Zorg dat we een lijst terugkrijgen
  )

  # Zet de lijst met losse stukjes terug om naar een geometry kolom (sfc)
  expanded_lines$geometry <- st_as_sfc(list_geoms, crs = st_crs(lines_sf))

  # ---------------------------------------------------

  # 6. Samenvoegen
  expanded_lines <- expanded_lines %>%
    select(-n_segments, -seg_num, -step_size, -start_frac, -end_frac)

  final_sf <- bind_rows(short_lines, expanded_lines) %>%
    arrange(orig_id)

  final_sf$orig_len <- NULL
  final_sf$orig_id <- NULL

  # Bereken nieuwe lengte
  final_sf$new_length <- as.numeric(st_length(final_sf))

  message(paste("Klaar! Resultaat:", nrow(final_sf), "segmenten."))
  return(final_sf)
}

fd_fine <- split_lines_equal_length(fd, max_length = 200)

build_river_network_auto_nodes <- function(lines_sf) {
  message("Graaf bouwen: data opschonen en knooppunten genereren...")

  # --- STAP 0: DATA WASSTRAAT (Cruciaal na knippen) ---

  # 1. Verwijder lege geometrieën
  lines_sf <- lines_sf[!st_is_empty(lines_sf), ]

  # 2. Als de data een mix is (GEOMETRYCOLLECTION), haal alleen de lijnen eruit
  # Dit repareert situaties waar punten en lijnen door elkaar staan
  if (any(grepl("GEOMETRY|COLLECTION", st_geometry_type(lines_sf)))) {
    lines_sf <- st_collection_extract(lines_sf, "LINESTRING")
  }

  # 3. Filter op geometry type: gooi alles weg dat nu nog steeds 'POINT' is
  # We willen enkel (MULTI)LINESTRINGS
  valid_types <- c("LINESTRING", "MULTILINESTRING")
  lines_sf <- lines_sf[st_geometry_type(lines_sf) %in% valid_types, ]

  # 4. Filter op lengte: alles met lengte 0 (of bijna 0) moet weg
  # Punten die per ongeluk lijnen zijn geworden (start=eind) geven problemen in de graaf
  lines_sf$check_len <- as.numeric(st_length(lines_sf))
  lines_sf <- lines_sf[lines_sf$check_len > 0.001, ] # Filter alles < 1mm
  lines_sf$check_len <- NULL # Opruimen

  # 5. Nu is het veilig om te casten naar LINESTRING
  lines_sf <- st_cast(lines_sf, "LINESTRING", warn = FALSE)

  message(paste(" -> Netwerk bevat", nrow(lines_sf), "geldige lijnsegmenten."))

  # ----------------------------------------------------

  # --- STAP 1: NODES BEPALEN ---

  # Start- en eindpunten bepalen
  p_start <- st_line_sample(lines_sf, sample = 0)
  p_end   <- st_line_sample(lines_sf, sample = 1)

  coords_start <- st_coordinates(p_start)
  coords_end   <- st_coordinates(p_end)

  # Maak Node ID's (afgerond op 1mm)
  node_id_start <- paste(round(coords_start[,1], 3), round(coords_start[,2], 3), sep = "_")
  node_id_end   <- paste(round(coords_end[,1], 3),   round(coords_end[,2], 3),   sep = "_")

  # --- STAP 2: ATTRIBUTEN ---

  lines_sf$from_node <- node_id_start
  lines_sf$to_node   <- node_id_end
  lines_sf$weight    <- as.numeric(st_length(lines_sf))

  # --- STAP 3: GRAAF BOUWEN ---

  edges <- st_drop_geometry(lines_sf) %>%
    select(from = from_node, to = to_node, weight, everything())

  g <- graph_from_data_frame(edges, directed = TRUE)

  message("Graaf succesvol gebouwd!")
  return(list(graph = g, network_sf = lines_sf))
}

river_network_fine <- build_river_network_auto_nodes(fd_fine)

head(beekprik)

#

get_beekprik_clusters <- function(observations_sf,
                                  network_list,
                                  target_years = NULL,
                                  max_dist_m = 5000,
                                  snap_tolerance = 25,
                                  time_col = "jaar") {

  message("1/5 Validatie en filteren van de data...")
  g <- network_list$graph
  net_sf <- network_list$network_sf

  # 1. Tijdsfilter toepassen (indien opgegeven)
  if (!is.null(target_years)) {
    obs_filtered <- observations_sf %>% filter(.data[[time_col]] %in% target_years)
  } else {
    obs_filtered <- observations_sf
  }

  if (nrow(obs_filtered) == 0) {
    warning("Geen observaties gevonden voor de opgegeven criteria.")
    return(NULL)
  }

  message(paste(" ->", nrow(obs_filtered), "observaties geselecteerd voor de analyse."))

  # 2. Snappen van observaties aan het netwerk
  message("2/5 Observaties snappen aan de graaf...")
  nearest_idx <- st_nearest_feature(obs_filtered, net_sf)
  dists <- st_distance(obs_filtered, net_sf[nearest_idx, ], by_element = TRUE)

  valid_snaps <- as.numeric(dists) <= snap_tolerance
  if (sum(valid_snaps) == 0) {
    warning("Geen enkele observatie lag binnen de snap_tolerance van het netwerk.")
    return(NULL)
  }

  # Haal de unieke startknooppunten op
  # We gebruiken de from_node van het dichtstbijzijnde segment
  snapped_nodes <- unique(as.character(net_sf$from_node[nearest_idx[valid_snaps]]))
  message(paste(" ->", length(snapped_nodes), "unieke netwerk-nodes geïdentificeerd als startpunten."))

  # 3. Bereken netwerkafstanden
  message(paste("3/5 Netwerk verkennen (max afstand:", max_dist_m, "m)..."))
  # mode = "all" is cruciaal: vissen kunnen stroomopwaarts én stroomafwaarts zwemmen
  # distances() geeft een matrix: Rijen = snapped_nodes, Kolommen = alle nodes in de graaf
  dist_matrix <- distances(g, v = snapped_nodes, mode = "all")

  # Voor elke node in de graaf zoeken we de minimale afstand tot EENDER WELKE observatie
  # apply met MARGIN=2 betekent: voer uit per kolom
  min_dists <- apply(dist_matrix, 2, min)

  # Selecteer alle nodes die binnen de drempelwaarde vallen
  valid_nodes <- colnames(dist_matrix)[is.finite(min_dists) & min_dists <= max_dist_m]

  # 4. Losse clusters identificeren via een sub-graaf
  message("4/5 Ruimtelijk gescheiden clusters identificeren...")
  # We maken een tijdelijke graaf die enkel bestaat uit de bereikbare nodes
  sub_g <- induced_subgraph(g, valid_nodes)

  # Bepaal de verbonden componenten (dit zijn onze losse areaalclusters!)
  comps <- components(sub_g)

  # Maak een tabelletje met Node ID en bijbehorend Cluster ID
  node_clusters <- data.frame(
    node_id = names(comps$membership),
    cluster_id = as.character(comps$membership)
  )

  # 5. Segmenten ophalen en samenvoegen tot MULTILINESTRING
  message("5/5 Geometrieën samenvoegen (Dissolve)...")

  # We koppelen de segmenten op basis van de from_node.
  # (Omdat elk segment in de sub-graaf verweven is, is dit voldoende)
  cluster_segments <- net_sf %>%
    inner_join(node_clusters, by = c("from_node" = "node_id"))

  # Groepeer per cluster en smelt de geometrieën samen
  final_clusters_sf <- cluster_segments %>%
    group_by(cluster_id) %>%
    summarise(
      n_segmenten = n(),
      totale_lengte_km = round(sum(weight) / 1000, 2), # weight was lengte in meters
      geometry = st_union(geometry), # Smelt lijnen samen
      .groups = "drop"
    ) %>%
    # Zorg dat de output netjes één of meerdere MULTILINESTRINGS zijn
    st_cast("MULTILINESTRING")

  message(paste("Klaar! Er zijn", nrow(final_clusters_sf), "aaneengesloten beekprik-clusters gegenereerd."))

  return(final_clusters_sf)
}

# Bereken de clusters (bijv. 3km bereik rondom observaties in recente jaren)

beekprik_clusters_05_14_3km <- get_beekprik_clusters(
  observations_sf = beekprik_aanwezig,
  network_list = river_network_fine,
  target_years = 2005:2014,   # Tijdsvenster
  max_dist_m = 3000,          # 3 kilometer netwerkafstand
  snap_tolerance = 25,        # Maximale foutmarge om punt op de lijn te plakken
  time_col = "jaar"
)

beekprik_2005_2014 <- beekprik_aanwezig %>% filter(jaar %in% 2005:2014)

beekprik_clusters_15_24_3km <- get_beekprik_clusters(
  observations_sf = beekprik_aanwezig,
  network_list = river_network_fine,
  target_years = 2015:2024,   # Tijdsvenster
  max_dist_m = 3000,          # 3 kilometer netwerkafstand
  snap_tolerance = 25,        # Maximale foutmarge om punt op de lijn te plakken
  time_col = "jaar"
)

beekprik_2015_2024 <- beekprik_aanwezig %>% filter(jaar %in% 2015:2024)


# Visualiseer met mapview
mapview(fd) +
  mapview(beekprik_clusters_15_24_3km, zcol = "cluster_id", lwd = 3) +
  mapview(beekprik_2015_2024, col.regions = "red")

mapview(fd) +
  mapview(beekprik_clusters_05_14_3km, color = "red", lwd = 3) +
  mapview(beekprik_clusters_15_24_3km, color = "green", lwd = 3)


###### methode met migratieknelpunten ######

# deze werkt nog niet naar behoren, of is alleszins super traag.

get_beekprik_clusters_tijdreeks <- function(observations_sf,
                                            network_list,
                                            knelpunten_sf = NULL,
                                            target_years = NULL,
                                            max_dist_m = 5000,
                                            snap_tolerance_obs = 50,
                                            snap_tolerance_knelpunt = 25,
                                            time_col = "jaar",
                                            vhas_col_knelpunt = "vhas",
                                            vhas_col_network = "VHAS") {

  # --- 1. Voorbereiding ---
  g_base <- network_list$graph
  net_sf <- network_list$network_sf

  if (!is.null(target_years)) {
    observations_sf <- observations_sf %>% filter(.data[[time_col]] %in% target_years)
  }

  if (nrow(observations_sf) == 0) {
    warning("Geen observaties gevonden voor de opgegeven criteria.")
    return(NULL)
  }

  unieke_jaren <- sort(unique(observations_sf[[time_col]]))
  message(paste("Analyse gestart voor", length(unieke_jaren), "verschillende jaren..."))

  # Lijst om resultaten per jaar in op te slaan
  all_clusters <- list()

  # --- 2. Iteratie per jaar ---
  for (yr in unieke_jaren) {
    message(paste0("\n--- Verwerken van jaar: ", yr, " ---"))

    # A. Data filteren voor huidig jaar
    obs_yr <- observations_sf %>% filter(.data[[time_col]] == yr)
    message(paste(" ->", nrow(obs_yr), "observaties."))

    # B. Observaties snappen
    idx_obs <- st_nearest_feature(obs_yr, net_sf)
    dist_obs <- as.numeric(st_distance(obs_yr, net_sf[idx_obs, ], by_element = TRUE))
    valid_obs <- dist_obs <= snap_tolerance_obs

    start_nodes <- unique(as.character(net_sf$from_node[idx_obs[valid_obs]]))
    if (length(start_nodes) == 0) {
      message(" -> Geen observaties succesvol gesnapt binnen de tolerantie. Overslaan.")
      next
    }

    # C. Knelpunten supersnel snappen (Correcte Logica)
    barrier_nodes <- character(0)

    if (!is.null(knelpunten_sf)) {
      knel_yr <- knelpunten_sf %>% filter(.data[[time_col]] == yr)

      if (nrow(knel_yr) > 0) {
        snapped_barriers <- rep(NA_character_, nrow(knel_yr))

        # Stap 1: Bereken de geometrische fallback voor ALLES in één keer (dit is het zware rekenwerk)
        fallback_idx <- st_nearest_feature(knel_yr, net_sf)
        fallback_dist <- as.numeric(st_distance(knel_yr, net_sf[fallback_idx, ], by_element = TRUE))
        valid_fallback <- fallback_dist <= 25 # Tolerantie van max 25m voor de fallback

        # Stap 2: Pas jouw specifieke regels toe per knelpunt
        for (i in seq_len(nrow(knel_yr))) {
          k_vhas <- knel_yr[[vhas_col_knelpunt]][i]
          match_found_via_vhas <- FALSE

          # Regel 1: Probeer te matchen op VHAS
          if (!is.na(k_vhas) && k_vhas != "") {
            net_idx <- which(net_sf[[vhas_col_network]] == k_vhas)

            if (length(net_idx) > 0) {
              # Zoek dichtstbijzijnde segment BINNEN de correcte VHAS
              rel_idx <- st_nearest_feature(knel_yr[i, ], net_sf[net_idx, ])
              abs_idx <- net_idx[rel_idx]

              # Snap hierop, ONGEACHT de afstand
              snapped_barriers[i] <- as.character(net_sf$from_node[abs_idx])
              match_found_via_vhas <- TRUE
            }
          }

          # Regel 2: Geen VHAS match gevonden? Gebruik de fallback (mét 25m restrictie)
          if (!match_found_via_vhas) {
            if (valid_fallback[i]) {
              snapped_barriers[i] <- as.character(net_sf$from_node[fallback_idx[i]])
            } else {
              # Ligt verder dan 25m en geen VHAS? Dan vervalt dit knelpunt
              snapped_barriers[i] <- NA
            }
          }
        }

        barrier_nodes <- unique(na.omit(snapped_barriers))
        message(paste(" ->", length(barrier_nodes), "knelpunten succesvol gesnapt."))
      } else {
        message(" -> Geen actieve migratieknelpunten in dit jaar.")
      }
    }

    # D. Graaf aanpassen voor het huidige jaar
    # Let op: als een observatie *exact* op een knelpunt valt, halen we die start_node eruit
    # anders faalt het algoritme. Hij verspreidt dan niet.
    start_nodes <- setdiff(start_nodes, barrier_nodes)

    if (length(start_nodes) == 0) {
      message(" -> Alle startpunten vallen exact op een knelpunt en zijn geblokkeerd. Overslaan.")
      next
    }

    # Verwijder de knooppunten uit de graaf. Hierdoor is de verbinding daar definitief doorgeknipt.
    g_yr <- g_base
    if (length(barrier_nodes) > 0) {
      # igraph filtert op vertex names
      nodes_to_delete <- V(g_yr)[name %in% barrier_nodes]
      g_yr <- delete_vertices(g_yr, nodes_to_delete)
    }

    # E. Netwerkafstanden en clusters berekenen
    dist_matrix <- distances(g_yr, v = start_nodes, mode = "all")
    min_dists <- apply(dist_matrix, 2, min)
    valid_nodes <- colnames(dist_matrix)[is.finite(min_dists) & min_dists <= max_dist_m]

    # Componenten bepalen op basis van bereikbare gebieden
    sub_g <- induced_subgraph(g_yr, valid_nodes)
    comps <- components(sub_g)

    node_clusters <- data.frame(
      node_id = names(comps$membership),
      cluster_id_temp = as.character(comps$membership)
    )

    # F. Samenvoegen tot MULTILINESTRING en jaar toevoegen
    cluster_segments <- net_sf %>%
      inner_join(node_clusters, by = c("from_node" = "node_id"))

    clusters_sf_yr <- cluster_segments %>%
      group_by(cluster_id_temp) %>%
      summarise(
        jaar = yr,
        n_segmenten = n(),
        totale_lengte_km = round(sum(weight) / 1000, 2),
        geometry = st_union(geometry),
        .groups = "drop"
      ) %>%
      mutate(cluster_id = paste0("Y", yr, "_C", cluster_id_temp)) %>%
      select(cluster_id, jaar, n_segmenten, totale_lengte_km, geometry) %>%
      st_cast("MULTILINESTRING")

    all_clusters[[as.character(yr)]] <- clusters_sf_yr
    message(paste(" ->", nrow(clusters_sf_yr), "clusters gegenereerd voor", yr))
  }

  # --- 3. Resultaten bundelen ---
  if (length(all_clusters) > 0) {
    final_clusters_sf <- bind_rows(all_clusters)
    message(paste("\nKlaar! In totaal", nrow(final_clusters_sf), "clusters gegenereerd over alle jaren heen."))
    return(final_clusters_sf)
  } else {
    message("\nKlaar! Geen clusters kunnen genereren.")
    return(NULL)
  }
}

beekprik_habitat_historiek <- get_beekprik_clusters_tijdreeks(
  observations_sf = beekprik,
  network_list = river_network_fine,
  knelpunten_sf = knelpunten_tijdreeks,  # Je nieuwe knelpuntenlaag
  target_years = 2015:2024,
  max_dist_m = 3000,
  snap_tolerance_obs = 50,
  snap_tolerance_knelpunt = 25,          # Strengere snap voor knelpunten
  time_col = "jaar",
  vhas_col_knelpunt = "vhas",            # Kolomnaam in je tibble
  vhas_col_network = "VHAS"              # Kolomnaam in fd
)

# Kijken naar de evolutie:
mapview(beekprik_habitat_historiek, zcol = "jaar", burst = TRUE)


get_beekprik_clusters_tijdreeks_fast <- function(observations_sf,
                                                 network_list,
                                                 knelpunten_sf = NULL,
                                                 target_years = NULL,
                                                 max_dist_m = 5000,
                                                 snap_tolerance_obs = 50,
                                                 time_col = "jaar",
                                                 vhas_col_knelpunt = "vhas",
                                                 vhas_col_network = "VHAS") {

  # --- 1. Voorbereiding ---
  g_base <- network_list$graph
  net_sf <- network_list$network_sf

  if (!is.null(target_years)) {
    observations_sf <- observations_sf %>% filter(.data[[time_col]] %in% target_years)
  }

  if (nrow(observations_sf) == 0) return(NULL)

  unieke_jaren <- sort(unique(observations_sf[[time_col]]))
  message(sprintf("Analyse gestart voor %d jaren...", length(unieke_jaren)))

  all_clusters <- list()

  # --- 2. Iteratie per jaar ---
  for (yr in unieke_jaren) {
    message(paste0("\n--- Verwerken van jaar: ", yr, " ---"))

    # A. Data filteren voor huidig jaar
    obs_yr <- observations_sf %>% filter(.data[[time_col]] == yr)

    # B. Observaties snappen
    idx_obs <- st_nearest_feature(obs_yr, net_sf)
    dist_obs <- as.numeric(st_distance(obs_yr, net_sf[idx_obs, ], by_element = TRUE))
    valid_obs <- dist_obs <= snap_tolerance_obs

    start_nodes <- unique(as.character(net_sf$from_node[idx_obs[valid_obs]]))
    if (length(start_nodes) == 0) {
      message(" -> Geen observaties gesnapt binnen tolerantie. Overslaan.")
      next
    }

    # C. Knelpunten supersnel snappen (Vectorised)
    barrier_nodes <- character(0)
    if (!is.null(knelpunten_sf)) {
      knel_yr <- knelpunten_sf %>% filter(.data[[time_col]] == yr)

      if (nrow(knel_yr) > 0) {
        # Stap 1: Doe de geometrische fallback voor ALLES tegelijk in 1 bewerking (razendsnel)
        fallback_idx <- st_nearest_feature(knel_yr, net_sf)
        snapped_barriers <- as.character(net_sf$from_node[fallback_idx])

        # Stap 2: Overschrijf enkel degene met een geldige VHAS-match
        # (Dit is snel omdat we enkel zoeken binnen een piepkleine subset van het netwerk)
        for (i in seq_len(nrow(knel_yr))) {
          k_vhas <- knel_yr[[vhas_col_knelpunt]][i]
          if (!is.na(k_vhas) && k_vhas != "") {
            net_idx <- which(net_sf[[vhas_col_network]] == k_vhas)
            if (length(net_idx) > 0) {
              rel_idx <- st_nearest_feature(knel_yr[i, ], net_sf[net_idx, ])
              snapped_barriers[i] <- as.character(net_sf$from_node[net_idx[rel_idx]])
            }
          }
        }
        barrier_nodes <- unique(na.omit(snapped_barriers))
        message(paste(" ->", length(barrier_nodes), "knelpunten gesnapt."))
      }
    }

    # D. Graaf knippen
    start_nodes <- setdiff(start_nodes, barrier_nodes)
    if (length(start_nodes) == 0) {
      message(" -> Alle startpunten geblokkeerd. Overslaan.")
      next
    }

    g_yr <- g_base
    if (length(barrier_nodes) > 0) {
      nodes_to_delete <- V(g_yr)[name %in% barrier_nodes]
      g_yr <- delete_vertices(g_yr, nodes_to_delete)
    }

    # E. SUPER NODE TRICK voor snelle routeberekening
    # Voeg een dummy node toe en verbind hem met afstand 0 aan alle start_nodes
    g_yr <- add_vertices(g_yr, 1, name = "SUPER_START")
    edges_to_add <- as.vector(rbind("SUPER_START", start_nodes))
    g_yr <- add_edges(g_yr, edges_to_add, weight = 0)

    # Reken nu 1 keer de afstanden uit vanaf de SUPER_START. Dit elimineert de matrix-berekening!
    dist_vector <- distances(g_yr, v = "SUPER_START", mode = "all")[1, ]

    # Filter de nodes die binnen de afstand vallen (en gooi de SUPER_START er weer uit)
    valid_nodes <- names(dist_vector)[dist_vector <= max_dist_m & names(dist_vector) != "SUPER_START"]

    # F. Samenvoegen tot MULTILINESTRING
    message(paste(" -> Geometry samenvoegen (bereik:", length(valid_nodes), "nodes)..."))

    sub_g <- induced_subgraph(g_yr, valid_nodes)
    comps <- components(sub_g)

    node_clusters <- data.frame(
      node_id = names(comps$membership),
      cluster_id_temp = as.character(comps$membership)
    )

    cluster_segments <- net_sf %>%
      inner_join(node_clusters, by = c("from_node" = "node_id"))

    clusters_sf_yr <- cluster_segments %>%
      group_by(cluster_id_temp) %>%
      summarise(
        jaar = yr,
        n_segmenten = n(),
        totale_lengte_km = round(sum(weight) / 1000, 2),
        geometry = st_union(geometry),
        .groups = "drop"
      ) %>%
      mutate(cluster_id = paste0("Y", yr, "_C", cluster_id_temp)) %>%
      select(cluster_id, jaar, n_segmenten, totale_lengte_km, geometry) %>%
      st_cast("MULTILINESTRING")

    all_clusters[[as.character(yr)]] <- clusters_sf_yr
  }

  # --- 3. Resultaten bundelen ---
  if (length(all_clusters) > 0) {
    final_clusters_sf <- bind_rows(all_clusters)
    message(paste("\nKlaar! In totaal", nrow(final_clusters_sf), "clusters gegenereerd."))
    return(final_clusters_sf)
  } else {
    return(NULL)
  }
}

beekprik_habitat_historiek <- get_beekprik_clusters_tijdreeks_fast(
  observations_sf = beekprik,
  network_list = river_network_fine,
  knelpunten_sf = knelpunten_tijdreeks,  # Je nieuwe knelpuntenlaag
  target_years = 2015,
  max_dist_m = 3000,
  snap_tolerance_obs = 50,
  # snap_tolerance_knelpunt = 25,          # Strengere snap voor knelpunten
  time_col = "jaar",
  vhas_col_knelpunt = "vhas",            # Kolomnaam in je tibble
  vhas_col_network = "VHAS"              # Kolomnaam in fd
)
