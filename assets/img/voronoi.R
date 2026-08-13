# Load required libraries
library(ggplot2)
library(sf)
library(deldir)
library(dplyr)
library(mgcv)

set.seed(1729)

# Global spatial bounds
xlim <- c(-10, 10)
ylim <- c(-10, 10)

# --- 1. Fit Surface from 50 Sampled Observations ---
N_obs <- 50
obs_x <- runif(N_obs, xlim[1], xlim[2])
obs_y <- runif(N_obs, ylim[1], ylim[2])

# DIAGONAL CONTINUUM + LOCAL SPATIAL BLOBS
base_wave <- 12 * sin(0.22 * (obs_x + obs_y))
blob_1    <-  9 * exp(-((obs_x - 3.5)^2 + (obs_y + 5)^2) / 8)  # Warm blob (bottom-right)
blob_2    <- -9 * exp(-((obs_x + 2.5)^2 + (obs_y - 1.0)^2) / 10) # Cool blob (top-left)
blob_3    <-  7 * exp(-((obs_x + 5.5)^2 + (obs_y + 0.5)^2) / 5)  # Blob straddling the diagonal

obs_z <- base_wave + blob_1 + blob_2 + blob_3 + rnorm(N_obs, mean = 0, sd = 7.5)
df_obs <- data.frame(x = obs_x, y = obs_y, z = obs_z)
# Filter observed points to top-left triangle only
df_obs_topleft <- df_obs %>% filter(y >= x)

# Fit continuous GAM surface model (k=30 allows enough flexibility for the blobs)
fit_gam <- gam(z ~ s(x, y, k = 30), data = df_obs)

# Predict smooth latent surface across a fine grid
grid_res <- 300
grid_df <- expand.grid(
  x = seq(xlim[1], xlim[2], length.out = grid_res),
  y = seq(ylim[1], ylim[2], length.out = grid_res)
)
grid_df$z_hat <- predict(fit_gam, newdata = grid_df)

# --- 2. Voronoi Surface & Organic Tile Partitioning ---
# Random seeds for the Voronoi geometry
N_voronoi <- 90 
vor_x <- runif(N_voronoi, xlim[1], xlim[2])
vor_y <- runif(N_voronoi, ylim[1], ylim[2])
df_vor_seeds <- data.frame(x = vor_x, y = vor_y)

# Map fine grid points to their nearest Voronoi seed
dists <- outer(grid_df$x, df_vor_seeds$x, `-`)^2 + outer(grid_df$y, df_vor_seeds$y, `-`)^2
grid_df$tile_id <- apply(dists, 1, which.min)

# Calculate average GAM-predicted value per Voronoi tile
tile_averages <- grid_df %>%
  group_by(tile_id) %>%
  summarise(z_avg = mean(z_hat), .groups = "drop")

# Build Voronoi polygons
vt <- deldir(df_vor_seeds$x, df_vor_seeds$y, rw = c(xlim[1], xlim[2], ylim[1], ylim[2]))
tiles <- tile.list(vt)

poly_list <- lapply(seq_along(tiles), function(i) {
  t_i <- tiles[[i]]
  st_polygon(list(cbind(c(t_i$x, t_i$x[1]), c(t_i$y, t_i$y[1]))))
})

voronoi_sf <- st_sfc(poly_list) %>%
  st_sf(tile_id = 1:length(tiles)) %>%
  left_join(tile_averages, by = "tile_id")

# Identify whole Voronoi tiles belonging to Bottom-Right (vor_y <= vor_x)
br_tile_ids <- which(vor_y <= vor_x)

# Select whole Voronoi tiles (organic mosaic edge)
voronoi_organic <- voronoi_sf %>% filter(tile_id %in% br_tile_ids)

# Top-Left Surface: Keep grid points belonging strictly to non-Bottom-Right tiles
grid_topleft <- grid_df %>% filter(!(tile_id %in% br_tile_ids))

# Extract the organic boundary line separating the two regions
organic_boundary <- st_boundary(st_union(voronoi_organic))

# --- 3. Plotting ---
library(scico)
fill_limits <- range(c(grid_topleft$z_hat, voronoi_organic$z_avg), na.rm = TRUE)

ggplot() +
  # Top-Left: Continuous GAM Predicted Surface
  geom_raster(data = grid_topleft, aes(x = x, y = y, fill = z_hat)) +
  
  # Top-Left: Observed Points Overlay
  geom_point(data = df_obs_topleft, aes(x = x, y = y),
             shape = 21, fill = NA, color = "black", 
             size = 2, stroke = 0.6) +
  
  # Bottom-Right: Whole Voronoi Polygons with Regional Means
  geom_sf(data = voronoi_organic, aes(fill = z_avg), 
          color = "black", linewidth = 0.25, inherit.aes = FALSE) +
  
  # Bold Organic Boundary Line replacing straight diagonal
  geom_sf(data = organic_boundary, color = "black", linewidth = 0.1, fill = NA) +
  
  # Shared Color Palette
  scale_fill_distiller(palette = "RdYlBu", limits = fill_limits) +
  
  # Square coordinates
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  
  # Clean square boundary framing without axes
  theme_void() +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.margin = margin(1, 1, 1, 1)
  )
ggsave("square_voronoi_surface.png", width = 4, height = 4, units = "in", dpi = 300, bg = "white")
