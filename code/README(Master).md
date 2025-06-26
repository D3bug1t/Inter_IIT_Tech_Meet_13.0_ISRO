# High Resolution Elemental Mapping of Lunar Surface

## Overview
The project leverages the CLASS dataset's SOLARANG header property,which specifies the angle between the sun′s position and the horizontal plane of the lunar surface from the orbiter's position. The concept of orbits arises from the observed continuous and periodic behavior of the solar angles, where it goes from a maxima to a particular minimum and back to the same maximum, marking the completion of an orbit. The critical angle (here named as back_angle) is taken as 91 degrees, above which the data is considered for background calculation for the flare files (files having SOLARANG < 90 degree) of that particular orbit. For every orbit the flare and background file names are stored in a csv file to use it ahead in calculating the initial guesses for gaussian fitting of the CLASS spectrum.

It incorporates the use of dataframe obtained from `extract_orbits` function, loads the flare and background files corresponding to each orbit and finds initial guesses (amplitude, center, sigma) of each element ('Mg', 'Al', 'Si', 'Fe', 'Ca', 'Ti', 'Cr', 'O').

Though the Galactic Cosmic-Ray (GCR) background is subtracted from the XRF data, it still needs further correction due to the scattered solar spectrum disturbing the actual spectral distribution following which the scattered solar spectrum is taken into account by fitting the background subtracted CLASS spectrum through power law and then normalizing the total background subtracted spectrum (Galactic Cosmic Ray and scattered solar spectrum components) with the power law fitted component. 

The `Spatial_Resolution.py` script performs spatial resolution mapping and visualization for the ratio of a given element to **Silicon (Si)**. It calculates grid-based averages for the specified ratio, using latitude and longitude bounding boxes, and visualizes the distribution over a lunar surface map.

---

## Requirements

- Python 3.x
- Pandas
- NumPy
- os 
- Astropy
- datetime
- scipy
- tqdm
- pyspeckit
- csv
- warnings
- Matplotlib (for visualizations)
- geopandas
- pillow
- dash-leaflet (for website)
- dash-extensions (for website)
- Access to year-wise folders having month wise CSV Files containing XRF FITS files data i.e. Longitudes, Latitudes and elemental ratios.

To install the required dependencies, use the following:

```bash
   pip install -r requirements.txt
```
## Running the DERM(Dynamic Elemental Ratio Mapping) software

In order to run the main DERM(Dynamic Elemental Ratio Mapping) terminal software, run the main.py in the terminal. 

```bash
   python3 main.py
```
Provide the path to the dataset to be used to the software(Only in the case when new data is to be processed, the data for the last 5 year has already been processed). If only the Analysis of the current data is needed you may just leave it blank.


## Running the Interactive Website

Run the website.py file with python and enter the desired element. This will host the website on local server. Open the link displayed to enter the website.

```bash
   python3 website.py
```

## Functions Overview

### `extract_orbits`

**Purpose:**  
Extracts the orbit using the CLASS files' solar angle (SOLARANG) property from its header, segregating the flare (having solar angle < 90 degree, thus considering them as daytime files) and background (having solar angle > 91 degree, thus considering them as background) files and storing the same in a .csv file.

**Parameters:**
- `data_path` (`str`): Path to the directory where CLASS files are stored.
- `back_angle` (`int`): Value of the solar angle above which the files are grouped as background files for that particular orbit,  taken as 91.
- `year` (`int`): Year of the CLASS files taken into consideration.
- `month` (`int`): Month of that particular year files taken into consideration.
- `day` (`int`): Day of that particlaur month files taken into consideration.

**Returns:**
- Returns a dataframe where each row depicts a particular orbit inside which the first column stores the flare files and the second column stores the background files of that orbit.

**Working:**  
1. The code appends a '0' infront of the date/month in case it is less than 10 (Example: 6 -> 06).
2. It loads the file of that particular date and observes the periodic behaviour of the solar angle to mark the completion of one orbit.
3. When it reaches the last peak of that day, it checks whether the path to the next day exists or not and accordingly appends that orbit to the dataframe.

---

### `find_nearest_peak`

**Purpose:**  
Finds the nearest spectral peak to a target energy level for a specified element.

**Parameters:**
- `final_peaks` (`array`): Array of peak intensities (amplitudes) in the spectrum.
- `final_channels` (`array`): Array of channel values (energies) corresponding to `final_peaks`.
- `element_channels` (`dict`): Dictionary mapping element names to their target energy levels in keV.
- `element` (`str`): Name of the element for which the peak is being searched (e.g., 'Mg', 'Al', 'Si', 'Fe', 'Ca', 'Ti', 'Cr', 'O').

**Returns:**
- `tuple`: Peak intensity and energy of the nearest peak if found, otherwise `(None, None)`.

**Working:**  
1. Looks up the target energy for the specified element from `element_channels`.
2. Computes the absolute difference between this target energy and each peak channel in `final_channels`.
3. Finds the nearest peak within a threshold (0.125 keV) of the target energy, returning the peak's intensity and channel.

---

### `coverage_catalogue`

**Purpose:**  
Processes XRF FITS files to analyze spectral data for a particular orbit at a time, identifies background subtracted spectral peaks for elements of interest, and outputs the results to a CSV.

**Parameters:**
- `abs_path` (`str`): The absolute path to the directory containing CLASS FITS files.
- `year` (`int`): Year of the CLASS files taken into consideration.
- `month` (`int`): Month of that particular year files taken into consideration.
- `output_file` (`str`): Path to the CSV file where the results will be saved.
- `significance_level` (`float`): Threshold value to determine significant peaks for elements, used to filter out noise and retain only peaks with sufficient prominence, here taken as 3.

**Working:**  
1. Segregates the flare and background files for a single orbit by iterating in the dataframe obtained from the extract_orbits function (based on the SOLARANG provided in the header of the CLASS file).
2. Background files:
   - Loads and averages the background counts from the background FITS files for that particular orbit.
3. Flare files:
   - Loads spectral data, confines it in an energy range of 0.5keV to 7keV and uses the prominence parameter of the find_peaks function to ensure dynamic peak handling of the CLASS spectrum.
4. Finds nearest peaks for target elements using `find_nearest_peak`.
5. Subtracts background and records amplitude and channel information for each significant peak, based on the `significance_level`, and saves the analysis results to the specified CSV file.

---

### `find_sigma`

**Purpose:**  
The widths of the peaks in the final spectrum are well approximated with Full Width at Half-Maximum (FWHM) which is used to calculate the sigma parameter of the initial guess while fitting the gaussian. 

**Parameters:**
- `filtered_channels` (array): Array of channels filtered in the range 0.5keV<=E<=7keV
- `filtered_counts` (array): Array of counts corresponding to `filtered_channels`.
- `el_amp` (float): Amplitude corresponding to the particular element from the database obtained from `initial_guesses`.
- `el_channel` (float): Channel of an element corresponding to `el_amp`. 

**Returns:**
- Returns the sigma parameter of the initial guesses while fitting the gaussian.

**Working:**  
1. Segregates the channels at half the amplitude of the guessed peak and finding the rightmost and leftmost channels corresponding to the same.
2. Using these obtained channels FWHM calculated is used to obtain the sigma parameter.

---

### `Gauss_fit`

**Purpose:**  
The final guesses generated(after removing the GCR background) are passed into this function to perform the spectral fitting.

**Parameters:**
- `filtered_channels` (array): Array of channels filtered in the range 0.5keV<=E<=7keV
- `filtered_counts` (array): Array of counts corresponding to `filtered_channels`.
- `guesses` (array): Array of initial guesses consisting of amplitude,center,sigma(in-order).

**Returns:**
- Returns the fitted parameters and uncertainities corresponding to the final values.

**Working:**  
It takes in the guessed parameters (amplitude, center, sigma) obtained from  the `initial_guesses` and `find_sigma` function for gaussian fitting after subtraction of both background components (Galactic cosmic ray and scattered solar spectrum) from the CLASS spectrum.

---

### `fitted_catalogue`

**Purpose:**  
Processes XRF FITS files to analyze spectral data for a particular orbit at a time, subtracts the galactic cosmic ray background component from the spectrum, incorporates power law fitting for taking into account the scattered solar spectrum and finally performs the multi-gaussian fit and outputs the results to a CSV.

**Parameters:**
- `abs_path` (`str`): The absolute path to the directory containing CLASS FITS files.
- `year` (`int`): Year of the CLASS files taken into consideration.
- `month` (`int`): Month of that particular year files taken into consideration.
- `output_file` (`str`): Path to the CSV file where the results will be saved.

**Returns:**
Returns the element-wise fitted parameters after multi-gaussian fitting as well as the scattered solar spectrum parameters.

**Working:**  
1. Segregates the flare and background files for a single orbit by iterating in the dataframe obtained from the extract_orbits function (based on the SOLARANG provided in the header of the CLASS file).
2. Background files: Loads and averages the background counts from the background FITS files for that particular orbit.
3. Flare files: Checks the occurrence of the start time in the prominence_file's `'Timestamp'` column and loads the spectral data.
4. Subtracts entire galactic cosmic ray background from the spectrum and checks the occurrence of not null peak amplitude values of the elements to compute the initial_guesses, and saves the analysis results to the specified CSV file.
5. Incorporates the scattered solar spectrum background component by power law fitting and subtracts the same from the spectrum obtained before.
6. Finally the spectrum is multi-gaussian fitted and the resutes are stored in a csv file.

---

### `normalize_peaks`

**Purpose:**
To normalizing the total background subtracted CLASS spectrum with respect to the power law background component.

**Parameters:**
- `path` (`str`): The absolute path to the directory containing multi-gaussian fitted csv files.
- `year` (`int`): Year of the CLASS files taken into consideration.
- `month` (`int`): Month of that particular year files taken into consideration.
- `output_file` (`str`): Path to the CSV file where the results will be saved.

**Returns:**
Returns the normalized element-wise parameters.

**Working:**
1. The multi-gaussian fitted catalog obtained from `fitted_catalog` function is loaded for extracting fitted peak amplitudes of the elements and the power law fitted scattered solar spectrum component.
2. The power law component normalizes the fitted peak amplitudes element-wise and saves the result in a csv file.

---

### `avg_km_per_degree_lon`

**Purpose:**
The function avg_km_per_degree_lon calculates the average distance (in kilometers) covered per degree of longitude for a specified range of latitudes on a spherical body (the Moon), given its radius.

**Parameters:**
- `lat1 (float)`: The starting latitude of the range (in degrees).
- `lat2 (float)`: The ending latitude of the range (in degrees).

**Returns:**
The function returns a `float`: The average distance (in kilometers) covered per degree of longitude for the given latitude range.

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

 ## Notebooks Overview:

 ### random_checking.ipynb
Check and validate the fitted spectrums

### Uncertainity.ipynb
Codes to generate uncertainity map used in report


---

