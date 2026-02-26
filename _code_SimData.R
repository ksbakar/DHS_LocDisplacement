
###############################################################################
###############################################################################

# R code to simulate data from the model 

library(sf)
library(dplyr)
library(ggplot2)
library(spatstat.random)

set.seed(3633)

################################################################################
# load - this code is for irregular geospatial areal shape 
# load ADM3 shapefile 
################################################################################

grid_sf <- read_sf("*ADM3.shp")
grid_sf$area_id <- 1:nrow(grid_sf)

# reproject polygon to UTM (meters) for sampling/displacement

grid_sf_proj <- st_transform(grid_sf, 32646)  # UTM zone 46N

# Sample original points inside polygon

n_total <- 35
n_urban <- 10
n_rural <- n_total - n_urban

# convert polygon to spatstat owin

grid_owin <- as.owin(st_geometry(st_union(grid_sf_proj)))

# sample n_total points

points_ppp <- spatstat.random::runifpoint(n_total, win = grid_owin)
coords <- data.frame(X = points_ppp$x, Y = points_ppp$y)

# sf in projected CRS

original_sf_proj <- st_as_sf(coords, coords = c("X","Y"), crs = 32646)
original_sf_proj$type <- c(rep("Urban", n_urban), rep("Rural", n_rural))

# back to WGS84

original_sf <- st_transform(original_sf_proj, 4326)

# displacement function

displace_points_fast <- function(points_sf, polygon_sf, radius_min = 0, radius_max = 2000) {
  
  # projected CRS
  
  polygon_proj <- st_transform(polygon_sf, 32646)
  points_proj <- st_transform(points_sf, 32646)
  
  n <- nrow(points_proj)
  coords <- st_coordinates(points_proj)
  displaced_coords <- matrix(NA, nrow = n, ncol = 2)
  
  for(i in 1:n) {
    valid <- FALSE
    tries <- 0
    while(!valid & tries < 50) {
      tries <- tries + 1
      angle <- runif(1, 0, 2*pi)
      distance <- sqrt(runif(1, radius_min^2, radius_max^2))
      new_x <- coords[i,1] + distance * cos(angle)
      new_y <- coords[i,2] + distance * sin(angle)
      new_point <- st_sfc(st_point(c(new_x, new_y)), crs = 32646)
      
      if(st_within(new_point, st_union(polygon_proj), sparse = FALSE)) {
        displaced_coords[i,] <- c(new_x, new_y)
        valid <- TRUE
      }
    }
    
    # keep original if not valid
    
    if(!valid) displaced_coords[i,] <- coords[i,]
  }
  
  displaced_sf <- st_as_sf(data.frame(lon_new = displaced_coords[,1],
                                      lat_new = displaced_coords[,2],
                                      type = points_sf$type),
                           coords = c("lon_new","lat_new"), crs = 32646)
  
  displaced_sf <- st_transform(displaced_sf, 4326)
  return(displaced_sf)
}

# displacement separately for Urban / Rural

urban_sf <- original_sf %>% filter(type=="Urban")
rural_sf <- original_sf %>% filter(type=="Rural")

displaced_urban <- displace_points_fast(urban_sf, grid_sf, radius_min=0, radius_max=2000)
displaced_rural <- displace_points_fast(rural_sf, grid_sf, radius_min=3000, radius_max=10000)

# combine

displaced_sf <- rbind(displaced_urban, displaced_rural)

# assign area_id

original_sf <- st_join(original_sf, grid_sf, join = st_within)
displaced_sf <- st_join(displaced_sf, grid_sf, join = st_within)

# compute area centroids

grid_midpoints <- st_centroid(grid_sf)

p1 <- ggplot() +
  geom_sf(data = grid_sf, fill = NA, color = "black", linetype = "dashed") +
  geom_sf(data = original_sf, aes(color = type, shape = type), size = 2, alpha = 1) +
  geom_sf(data = displaced_sf, shape = 3, size = 3, alpha = 1) +
  geom_sf(data = grid_midpoints, shape = 7, color = "black", size = 2) +
  labs(color = "", shape="") +
  theme_classic() +
  theme(legend.position = 'bottom')

################################################################################
# load - this code is for grid areal data 
################################################################################

library(sf)
library(ggplot2)
library(dplyr)
library(tidyr)

set.seed(3633)

# displace points (unchanged)

displace_points <- function(lon, lat, n, radius_min = 0, radius_max = 1000, bbox = NULL) {
  earth_radius <- 6371000
  displaced <- data.frame(lon_new = numeric(n), lat_new = numeric(n))
  
  for(i in 1:n) {
    repeat {
      angle <- runif(1, 0, 2*pi)
      distance <- sqrt(runif(1, radius_min^2, radius_max^2))
      
      delta_lat <- distance * cos(angle) / earth_radius
      delta_lon <- distance * sin(angle) / (earth_radius * cos(lat * pi / 180))
      
      new_lat <- lat + delta_lat * (180/pi)
      new_lon <- lon + delta_lon * (180/pi)
      
      if(is.null(bbox) || (new_lon >= bbox["xmin"] & new_lon <= bbox["xmax"] &
                           new_lat >= bbox["ymin"] & new_lat <= bbox["ymax"])) {
        displaced$lon_new[i] <- new_lon
        displaced$lat_new[i] <- new_lat
        break
      }
    }
  }
  return(displaced)
}

# original points

n_total <- 35
n_urban <- 10
n_rural <- n_total - n_urban

lat_center <- 40.0
lon_center <- -75.0
lat_range <- 0.1
lon_range <- 0.1

original_points <- data.frame(
  lon = runif(n_total, lon_center - lon_range, lon_center + lon_range),
  lat = runif(n_total, lat_center - lat_range, lat_center + lat_range),
  type = c(rep("Urban", n_urban), rep("Rural", n_rural))
)

original_sf <- st_as_sf(original_points, coords = c("lon", "lat"), crs = 4326)
bbox <- st_bbox(original_sf)

# displace points

displaced_points <- original_points %>%
  rowwise() %>%
  mutate(
    new_coords = list(
      if(type == "Urban") {
        displace_points(lon, lat, 1, 0, 2000, bbox)
      } else {
        displace_points(lon, lat, 1, 3000, 10000, bbox)
      }
    )
  ) %>%
  unnest(cols = c(new_coords))

displaced_sf <- st_as_sf(displaced_points, coords = c("lon_new", "lat_new"), crs = 4326)

# 12-area grid (3x4)

all_points <- rbind(
  st_geometry(original_sf) %>% st_sf(type = original_sf$type),
  st_geometry(displaced_sf) %>% st_sf(type = displaced_sf$type)
)
bbox <- st_bbox(all_points)

grid <- st_make_grid(st_as_sfc(bbox), n = c(3, 4), what = "polygons")
grid_sf <- st_sf(area_id = 1:length(grid), geometry = grid)

# area_id to points using st_within

original_sf <- st_join(original_sf, grid_sf, join = st_within)
displaced_sf <- st_join(displaced_sf, grid_sf, join = st_within)

# fix points that are NA by nearest area centroid

assign_nearest <- function(points_sf, grid_sf) {
  na_idx <- which(is.na(points_sf$area_id))
  if(length(na_idx) > 0) {
    centroids <- st_centroid(grid_sf)
    nearest_idx <- st_nearest_feature(points_sf[na_idx, ], centroids)
    points_sf$area_id[na_idx] <- centroids$area_id[nearest_idx]
  }
  return(points_sf)
}

original_sf <- assign_nearest(original_sf, grid_sf)
displaced_sf <- assign_nearest(displaced_sf, grid_sf)

# midpoints of each area

grid_midpoints <- st_centroid(grid_sf)

# plot

p2 <- ggplot() +
  geom_sf(data = grid_sf, fill = NA, color = "black", linetype = "dashed") +
  geom_sf(data = original_sf, aes(color = type, shape = type), size = 2, alpha = 1) +
  geom_sf(data = displaced_sf, shape = 3, size = 3, alpha = 1) +
  geom_sf(data = grid_midpoints, shape = 7, color = "black", size = 2) +
  labs(color = "", shape="") +
  theme_classic() +
  theme(legend.position = 'bottom') 

#print(seed)

gridExtra::grid.arrange(p1,p2,ncol=2)

# data simulation (code is similar for both grid and non-grid areal units)

library(sf)
library(dplyr)
library(Matrix)
library(geodist)

set.seed(123)

J   <- nrow(original_sf)   # 35 clusters
m   <- nrow(grid_sf)       # 12 areas
n_j <- 10
N   <- J * n_j

# coordinates

cluster_coords_true <- st_coordinates(original_sf)   # J x 2 (lon, lat)
cluster_coords_disp <- st_coordinates(displaced_sf)
area_coords <- st_coordinates(grid_midpoints)        # m x 2

area_id <- original_sf$area_id
urban <- ifelse(original_sf$type == "Urban", 1, 0)

# geodetic distances (km)

dist_true <- geodist(
  cluster_coords_true[, c("X","Y")],
  area_coords[, c("X","Y")],
  measure = "geodesic"
) / 1000   # meters -> km

# urban and rural displacement noise (km)

sigma_u <- 1 #0.5   # urban SD km
sigma_r <- 3 #2.0   # rural SD km

eps <- ifelse(
  urban == 1,
  rnorm(J, 0, sigma_u),
  rnorm(J, 0, sigma_r)
)

dist_disp <- dist_true + eps
dist_disp[dist_disp < 0] <- 0.001

# phi

phi_fun <- function(d, range = 50) { 
  exp(-d / range)
}

Phi_star <- phi_fun(dist_disp)   # J x m

# moran operator M(W)

W <- matrix(0, J, J)
for (k in unique(area_id)) {
  idx <- which(area_id == k)
  W[idx, idx] <- 1
}
diag(W) <- 0

Z <- matrix(1, J, 1)
Pz <- Z %*% solve(t(Z) %*% Z) %*% t(Z)

M <- (diag(J) - Pz) %*% W %*% (diag(J) - Pz)

# latent spatial process

# areal mean

mu_z_area <- rnorm(m, mean = 0, sd = 0.5)
mu_z <- mu_z_area[area_id]

# latent coefficients

eta_tilde <- rnorm(m)

Psi <- M %*% Phi_star
z_star <- as.vector(mu_z + Psi %*% eta_tilde)

# observed spatial covariate z_D

sigma_nu <- 0.1 #
z_D <- z_star + rnorm(J, 0, sigma_nu)

# individual-level

x_ij <- rnorm(N)

alpha <- 1
beta  <- 0.5
zeta  <- 0.5

linpred <- alpha +
  beta * x_ij +
  zeta * rep(z_star, each = n_j)

p <- plogis(linpred)
y <- rbinom(N, 1, prob = p)

dat_sim <- data.frame(
  cluster_id = rep(1:J, each = n_j),
  area_id    = rep(area_id, each = n_j),
  urban      = rep(urban, each = n_j),
  x_ij       = x_ij,
  z_star     = rep(z_star, each = n_j),
  z_D        = rep(z_D, each = n_j),
  y          = y
)

head(dat_sim)
summary(dat_sim)

###############################################################################
###############################################################################