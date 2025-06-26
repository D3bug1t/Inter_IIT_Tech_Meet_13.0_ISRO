#  Calculating the Initial Guesses for gaussian fitting of the CLASS spectrum

## Overview

This code incorporates the use of dataframe obtained from `extract_orbits` function, loads the flare and background files corresponding to each orbit and finds initial guesses (amplitude, center, sigma) of each element ('Mg', 'Al', 'Si', 'Fe', 'Ca', 'Ti', 'Cr', 'O').

## Functions Overview

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