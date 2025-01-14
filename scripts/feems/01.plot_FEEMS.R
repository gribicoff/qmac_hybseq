## Load packages
library(scales)
library(readr)
library(sf)
library(magrittr)
library(dplyr)
library(maps)

## Define functions
prepare_edge_data <- function(nodes_df, edges_df) {
  edges_df %>%
    mutate(
      from = as.numeric(from),
      to = as.numeric(to)
    ) %>% 
    # Join coordinates for node1
    left_join(nodes_df, 
              by = c("from" = "node_id")) %>%
    rename(node1_lat = latitude,
           node1_lon = longitude) %>%
    # Join coordinates for node2
    left_join(nodes_df,
              by = c("to" = "node_id")) %>%
    rename(node2_lat = latitude,
           node2_lon = longitude)
}

map_values_to_colors <- function(values, colors = c("red1", "dodgerblue"), n_colors = 100) {
  # Create color palette
  color_palette <- colorRampPalette(colors)(n_colors)
  
  # Normalize values to 0-1 range
  normalized_values <- (values - min(values)) / (max(values) - min(values))
  
  # Map normalized values to color indices
  color_indices <- floor(normalized_values * (n_colors - 1)) + 1
  
  # Return colors
  return(color_palette[color_indices])
}

## Read data
# U.S. Shapefile
us = st_read("OUT/state.shp")
# Edge locations
edges = read_csv(file = "OUT/edges_weights.csv") %>% 
  rename(from = node1,
         to = node2,
         edge_value = weight) %>% 
  mutate(from = as.character(from),
         to = as.character(to)) %>% 
  dplyr::rename()
# Node locations
nodes = read_csv(file = "OUT/node_positions.csv") %>% 
  mutate(ID = as.numeric(0:(nrow(.)-1))) %>% 
  select(ID, longitude, latitude) %>% 
  rename(node_id = ID)
# Sample locations
samples = read_csv(file = "OUT/sample_positions.csv") %>% 
  st_as_sf(coords = c("longitude", "latitude")) %>% 
  st_coordinates() %>% 
  as.data.frame()

## Format data
edges_with_coords = prepare_edge_data(nodes_df = nodes, 
                                      edges_df = edges)

# Set up edge colors and colors for legend
cols = map_values_to_colors(values = log(edges_with_coords$edge_value))
legend_cols = map_values_to_colors(values = log(edges_with_coords$edge_value)[order(log(edges_with_coords$edge_value))])

# Set up map limits
xlim = range(nodes$longitude)
ylim = range(nodes$latitude)
lims = range(log(edges_with_coords$edge_value))

# Formatting
lwd = 1 - log(edges_with_coords$edge_value)/max(log(edges_with_coords$edge_value))
alpha = 1 - log(edges_with_coords$edge_value)/max(log(edges_with_coords$edge_value))/max(1 - log(edges_with_coords$edge_value)/max(log(edges_with_coords$edge_value)))
col = alpha(colour = cols, 
            alpha = alpha)

# Plot map
# pdf(file = "OUT/20250113_QMAC_FEEMS.pdf",
#     width = 4.65,
#     height = 5.21)
plot(
  st_geometry(us),
  xlim = xlim,
  ylim = ylim,
  border = "grey40",
  lwd = 1
)
points(
  x = nodes$longitude,
  y = nodes$latitude,
  pch = 20,
  cex = 0.10,
  las = 1,
  col = "black")
segments(
  x0 = edges_with_coords$node1_lon,
  y0 = edges_with_coords$node1_lat,
  x1 = edges_with_coords$node2_lon,
  y1 = edges_with_coords$node2_lat,
  col = col,
  lwd = lwd
)
box(lwd = 2)
axis(side = 1)
axis(side = 2, las = 1)
points(
  x = samples$X,
  y = samples$Y,
  pch = 20,
  cex = 1
)
legend("topleft",
       legend = "Populations",
       pch = 20,
       cex = 1,
       bty = 'n')
phytools::add.color.bar(cols = legend_cols,
                        leg = 40,
                        title = expression(paste("Log"[10],"(w)")),
                        subtitle = "",
                        lims = range(log(edges_with_coords$edge_value)),
                        prompt = F,
                        x = -105,
                        y = 21.5,
                        outline = F,
                        lwd = 7.5)
# dev.off()
















