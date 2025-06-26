# Extracting Orbits for calculating appropriate initial guesses

## Overview

The project leverages the CLASS dataset's SOLARANG header property,which specifies the angle between the sun′s position and the horizontal plane of the lunar surface from the orbiter's position. The concept of orbits arises from the observed continuous and periodic behavior of the solar angles, where it goes from a maxima to a particular minimum and back to the same maximum, marking the completion of an orbit. The critical angle (here named as back_angle) is taken as 91 degrees, above which the data is considered for background calculation for the flare files (files having SOLARANG < 90 degree) of that particular orbit. For every orbit the flare and background file names are stored in a csv file to use it ahead in calculating the initial guesses for gaussian fitting of the CLASS spectrum.

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
