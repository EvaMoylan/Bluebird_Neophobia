################################################################################
# Step 2: LiDAR Processing for Bluebird Neophobia Nestboxes
# Updated: 10/25/26
# Objective: Calculate habitat structure at nest boxes from LiDAR data
# Authors: Lindsay Peria (template code), Isaiah Woodard (template code, .LAZ data), Eva Moylan (project-specific code)
################################################################################

# This won't work with current directory - LAZ files too large so stored on Mac
# Script will need to be updated when the missing box is added 

#################################### Set-up #################################### 

# Load packages #
library(tidyverse)
library(lidR)
library(terra)
library(sf)
library(here)

# Load .LAZ files containing compressed LiDAR point cloud data for the entire field system (from Isaiah Woodard)
las.cat <- readLAScatalog(here("Data", "LAZ_files")) # Can ignore overlap warning

class(las.cat) # Checks class
crs(las.cat) # Checks CRS

# Load points for all nest boxes
box_points <- st_read(here("Data","Nestbox_points","Nestbox_points.shp"))
box_points <- st_transform(box_points, crs = 26913) # Match CRS to LAZ files

# Load a key for boxes that received bluebird neophobia trials
neo <- read_csv(here("Data", "neo_boxes.csv")) |>
  rename(Name = nest_box) # Rename for easy joining with spatial dataset

# Subset box .shp file to include only necessary boxes
neo_points <- box_points |>
  semi_join(neo, by = "Name") |>
  mutate(Number = 1:length(unique(Name)),
         .before = Name)

# Looks like one was missing, need to add later
anti_join(neo, neo_points, by = "Name") # TP07B14

# Create an empty dataframe for the LiDAR output
lidar <- data.frame(matrix(NA, ncol = 4, nrow = 42)) # Change to 43 when TP07B14 is added

# Create empty columns in the LiDAR dataframe
names(lidar) <- c("point", "sd_height", "perchalfm", "rumple")

lidar$point <- 1:42

######################## Calculate Percent Canopy Cover ######################## 

# Canopy cover/ tree densities are what likely impact detection of first approach
# Calculate cover above 1.5 meters to remove shrubs and ground cover that might skew outputs

# Create a vector to catch failed calculations for nest boxes
failed_points <- integer(0)

for (i in neo_points$Number) { # For each nest box
  
  message("Processing nest box ", i) # Processing message in console
  
  tryCatch({
    
    # Clip area of interest
    roi <- clip_roi(las.cat, neo_points[i, ], radius = 40) # 50 meters is the average distance 
    
    # Guard against empty clips
    if (is.empty(roi)) {
      stop("Empty Region of Interest")
    }
    
    # Guard against missing ground points
    if (length(filter_ground(roi)) == 0) {
      stop("No ground points")
    }
    
    # Terrain and canopy models
    dtm <- rasterize_terrain(
      roi,
      res = 1,
      algorithm = tin()
    )
    
    dcm <- rasterize_canopy(
      roi,
      res = 1,
      algorithm = p2r(na.fill = tin())
    )
    
    # Rumple index
    lidar[i, "rumple"] <- rumple_index(dcm)
    
    # Canopy height
    height <- dcm - dtm
    
    # SD of canopy height
    lidar[i, "sd_height"] <- global(height, fun = "sd", na.rm = TRUE)
    
    # Percent canopy > 1.5 m
    tall <- freq(height >= 1.5)
    
    # Extract count for TRUE (value == 1) or default to 0 if no canopy exists
    canopy_pixels <- tall$count[tall$value == 1]
    if (length(canopy_pixels) == 0) canopy_pixels <- 0
    
    lidar[i, "perc_can"] <- (canopy_pixels / sum(tall$count)) * 100
    
  }, error = function(e) {
    
    message("⚠️  Failed at nest box ", i, ": ", e$message)
    
    failed_points <<- c(failed_points, i)
    
    # Leave row as NA and continue
    return(NULL)
  })
}

# Inspect failures
unique(failed_points) # No failed points

# If failed: subset points that failed
failed_nestboxes <- points_subset[failed_points, ]
lidar_nestbox <- lidar_empty |>
  left_join(points_subset,
            by = join_by("point" == "Number")) %>% 
  select(Name, point, sd_height, perc_can, rumple)

# Add nest box IDs back into the dataframe
final_lidar <- lidar |>
  left_join(neo_points[,c("Number","Name")], join_by("point" == "Number")) |>
  st_drop_geometry() |>
  select(-geometry) |>
  rename(nest_box = Name) |>
  relocate(nest_box, .after = point)

# Save file
write_csv(final_lidar, here("Outputs","lidar_metrics_40m_canopy.csv"))

