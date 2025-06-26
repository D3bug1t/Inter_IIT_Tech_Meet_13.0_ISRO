import os
import csv
from Extract_Orbits import *
from Initial_Guesses import *
from Power_law_fitting import *
from Spatial_Resolution import *
import datetime as dt
import warnings

warnings.filterwarnings("ignore")

def show_welcome_message():
    print("=======================================")
    print("           Welcome to DERM!            ")
    print("    Dynamic Elemental Ratio Mapping    ")
    print("=======================================")
    print("            Version 1.0.1              ")
    print("        Developed by Team 30           ")
    print("=======================================")


show_welcome_message()


abs_path = input("Provide The Path to directory listing CLASS fits files:- ") ## Example:- /path/to/your/CLASS/
years = list(input("Provide the years for Analysis(Comma-Separated):- ").split(','))
back_angle = 91
flag = 'y'
while(True and flag.lower()=='y'):
    flag1 = input("Proceed with Orbit Extraction(for background and flare categorisation)[Y/N]? ")
    flag2 = input("Proceed with Coverage Catalogue[Y/N]? ")
    flag3 = input("Proceed with Spectral Fitting[Y/N]? ")
    flag4 = input("Proceed with Ratio Calculation and Mapping[Y/N]? ")
    folders = os.listdir(os.getcwd())
    available_data = [folder for folder in folders if folder.isnumeric()]
    print(available_data)
    for year in years:
        if flag1.lower()=='y' and year not in available_data:
            for month in tqdm(range(1,13)):  # Loop through months 1 to 12
                month_str = f"{month:02d}"  # Zero-pad month as '01', '02', ..., '12'
                month_path = os.path.join(abs_path, f"ch2_cla_l1_{year}_{month_str}/cla/data/calibrated/", str(year), month_str)
                # Check if month directory exists
                if not os.path.exists(month_path):
                    print(f"Month directory not found: {month_path}")
                    continue

                # Loop through all days in the month
                for day in range(1, 32):  # Loop through days 1 to 31
                    try:
                        # Check if the date is valid
                        date = dt.datetime(int(year), int(month), int(day))
                        day_str = f"{day:02d}"  # Zero-pad day as '01', '02', ..., '31'
                        day_path = f"{month_path}/{day_str}"

                        # Check if the day's directory exists
                        if not os.path.exists(day_path):
                            print(f"Day directory not found: {day_path}")
                            continue

                        # Process the data for the day
                        orbits = extract_orbits(abs_path, back_angle, int(year), int(month), int(day))

                        # Save the extracted data to a CSV file
                        output_dir = f"{year}/{month_str}"
                        os.makedirs(output_dir, exist_ok=True)  # Create output directory if not exists
                        output_file = f"{output_dir}/{day_str}.csv"

                        with open(output_file, "w", newline="") as file:
                            writer = csv.writer(file)
                            writer.writerows(orbits)

                        print(f"Data processed and saved for {date.date()}")

                    except ValueError:
                        # Skip invalid dates (e.g., February 30, etc.)
                        print(f"Invalid date: {year}-{month:02d}-{day:02d}")
                        continue
        folders = os.listdir(os.path.join(os.getcwd(),'xrf_line_catalog'))
        if flag2.lower()=='y' and year not in folders:
            ## Coverage Catalog
            for month in tqdm(range(1,13)):  # Loop through months 1 to 12
                month_str = f"{month:02d}"  # Zero-pad month as '01', '02', ..., '12'
                month_path = f'{year}/{month_str}'
                if not os.path.exists(month_path):
                    print(f"Month directory not found: {month_path}")
                    continue
                coverage_catalogue(abs_path,year,month_str, f'xrf_line_catalog/{year}/{month_str}.csv',2)
            print("CSVs created for Initial Guesses at xrf_line_catalog/")

        if not os.path.exists(f'xrf_line_catalog_power_law_normalized'):
                        os.makedirs(f'xrf_line_catalog_power_law_normalized')

        folders = os.listdir(os.path.join(os.getcwd(),'xrf_line_catalog_power_law_normalized'))
        if flag3.lower()=='y' and year not in folders:
            ## Fitted Catalog
            for month in tqdm(range(1,13)):  # Loop through months 1 to 12
                month_str = f"{month:02d}"  # Zero-pad month as '01', '02', ..., '12'
                month_path = f'{year}/{month_str}'
                if not os.path.exists(month_path):
                    print(f"Month directory not found: {month_path}")
                    continue
                fitted_catalogue(abs_path,year,month_str, f'xrf_line_catalog_power_law/{year}/{month_str}.csv')
                if os.path.exists(f"xrf_line_catalog_power_law/{year}/{month_str}.csv"):
                    normalize_peaks(f'xrf_line_catalog_power_law/', year, month_str, output_file=f'xrf_line_catalog_power_law_normalized/{year}/{month_str}.csv' )
            print("CSVs created for Normalized Fitted Catalog with Background Removal at xrf_line_catalog_power_law_normalized/")
    if flag4.lower()=='y':
        element = input("Provide the element to be analysed:- ") ## Example:- Mg,Al,Si
        years_map = list(input("Provide the years for mapping the elemental ratios(Comma-Separated):- ").split(','))
        contours = input("Do the maps be plotted with Mare Contours [Y/N]? ")
        resolution = float(input("Specify the grid Resolution desired [in Km][Ex:-12.5,6.25,3.125 etc.]:- ")) 
        if(contours.lower()=='y'):
            Ratio_Mapping_with_contours(years_map,element,resolution)
        else:
            Ratio_Mapping_without_contours(years_map,element,resolution)
    flag5 = input("Proceed with creating shapefile[Y/N]? ")
    elements = list(input("Please provide the elements for which shapefiles are to be created(Comma-Separated):-").split(','))
    if not os.path.exists(f'shape_files'):
            os.makedirs(f'shape_files')
    if(flag5.lower()=='y'):
        for element in elements:
            print(elements)
            csv_path = os.path.join(os.getcwd(),"ratios",f"grid_average_{element}_values.csv")
            create_shp(csv_path,element)

    flag = input("Concat Another Year?[Y/N]?")
    if(flag.lower()=='y'):
        year = list(input("Provide the year for Analysis(Comma-Separated):- ").split(','))
    



