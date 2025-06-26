# Spatial Resolution Of Elements

## Overview
The `Spatial_Resolution.py` script performs spatial resolution mapping and visualization for the ratio of a given element to **Silicon (Si)**. The script calculates grid-based averages for the specified ratio, using latitude and longitude bounding boxes, and visualizes the distribution over a lunar surface map.


## Functions Overview

### avg_km_per_degree_lon(lat1,lat2)

**Purpose:**
The function avg_km_per_degree_lon calculates the average distance (in kilometers) covered per degree of longitude for a specified range of latitudes on a spherical body (the Moon), given its radius.

**Parameters:**
- `lat1 (float)`: The starting latitude of the range (in degrees).
- `lat2 (float)`: The ending latitude of the range (in degrees).

**Returns:**
The function returns a `float`:
- The average distance (in kilometers) covered per degree of longitude for the given latitude range.

**Working:** 
1. `Latitude to radians`:
 Convert the latitude values from degrees to radians using np.radians.
 The cosine of the latitude is used because the longitude lines converge at the poles, making the effective distance per degree of longitude proportional to the cosine of the latitude.
2. `Average cosine value`:
 The function calculates the cosine of the latitude for 100 evenly spaced points between lat1 and lat2 using np.linspace and np.mean.
3. `Distance calculation`:
 The circumference of the spherical body at the equator is calculated as 2π×radius.
 This circumference is divided by 360 to get the distance covered per degree of longitude at the equator.
 The value is multiplied by the average cosine of the latitudes to account for the variation in longitude distance with latitude.

---
### `Ratio_Mapping_with_contours`

**Purpose:**
 - Maps the ratio of a specified element to Si on a grid over the lunar surface.

**Parameters:**
- `year`: List of folders containing csv files which contain processed csv files month -wise.
- `element`: Element for which the ratio with Si is calculated.
- `resolution`: Desired resolution for the mapping to be done(in Km).

**Returns:**
Generates a visual representation of the grid on the lunar surface with Mare Boundaries.

**Working:**
1. `Data Collection`:
 - Iterates through the specified folders and reads all CSV files.
 - Combines the data into a single DataFrame.
 - Filters the data to ensure valid longitude values and non-zero   amplitudes for Si and the specified element.
2. `Ratio Calculation`:
 - Calculates the ratio of the specified element's peak amplitude to Si's peak amplitude for each data entry.
 - Filters out ratios exceeding 2 for better visualization and noise reduction
3. `Grid Creation`:
 - Divides the lunar surface into latitude and longitude bins based on the specified grid resolution.
 - Uses the avg_km_per_degree_lon function to calculate longitude bin sizes.
 - Assigns the calculated element-to-Si ratios to grid cells and calculates average ratio for each grid cell.4
4. `Visualisation`:
 - Generates a visual representation of the grid on the lunar surface.

---
### `Ratio_Mapping_without_contours`

**Purpose:**
 - Maps the ratio of a specified element to Si on a grid over the lunar surface.

**Parameters:**
- `year`: List of folders containing csv files which conatin processed csv files month -wise.
- `element`: Element for which the ratio with Si is calculated.

**Returns:**
Generates a visual representation of the grid on the lunar surface.

**Working:**
1. `Data Collection`:
 - Iterates through the specified folders and reads all CSV files.
 - Combines the data into a single DataFrame.
 - Filters the data to ensure valid longitude values and non-zero   amplitudes for Si and the specified element.
2. `Ratio Calculation`:
 - Calculates the ratio of the specified element's peak amplitude to Si's peak amplitude for each data entry.
 - Filters out ratios exceeding 2 for better visualization and noise reduction
3. `Grid Creation`:
 - Divides the lunar surface into latitude and longitude bins based on the specified grid resolution.
 - Uses the avg_km_per_degree_lon function to calculate longitude bin sizes.
 - Assigns the calculated element-to-Si ratios to grid cells and calculates average ratio for each grid cell.4
4. `Visualisation`:
 - Generates a visual representation of the grid on the lunar surface without Mare boundaries.

---