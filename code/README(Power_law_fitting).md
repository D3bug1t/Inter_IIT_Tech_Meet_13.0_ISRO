# Solar Continuum Modelling and Subtraction using Power Law

## Overview

Though the Galactic Cosmic-Ray (GCR) background is subtracted from the XRF data, it still needs further correction due to the scattered solar spectrum disturbing the actual spectral distribution. The following script contains the incorporation of this scattered solar spectrum by fitting the background subtracted CLASS spectrum through power law and then normalizing the total background subtracted spectrum (Galactic Cosmic Ray and scattered solar spectrum components) with the power law fitted component. 

## Functions Overview

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