from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os
from datetime import datetime, timedelta
import pandas as pd
import shutil
from scipy.special import wofz
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.stats import tstd
import csv
import shapefile
from tqdm import tqdm
from scipy.signal import find_peaks
import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon, Point
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
# import Create_shp.ipynb

def create_shp(csv_path,element):
        
    
    # Step 1: Read the CSV file
    csv_file = csv_path  # Replace with the path to your CSV file
    df = pd.read_csv(csv_file)

    # Step 2: Create geometries for the grid
    # Assuming the CSV contains `lat_min`, `lat_max`, `lon_min`, `lon_max` for grid boundaries
    geometries = []
    for _, row in df.iterrows():
        lat_min, lat_max = row["Latitude_Min"], row["Latitude_Max"]
        lon_min, lon_max = row["Longitude_Min"], row["Longitude_Max"]
        # Create a polygon for the grid cell
        geometries.append(
            Polygon([
                (lon_min, lat_min),
                (lon_min, lat_max),
                (lon_max, lat_max),
                (lon_max, lat_min),
                (lon_min, lat_min)  # Close the polygon
            ])
        )

    # Create a GeoDataFrame
    gdf = gpd.GeoDataFrame(df, geometry=geometries, crs="EPSG:4326")  # Use EPSG:4326 for lat/lon
    gdf["Grid_ID"] = gdf.apply(lambda row: f"Grid_{row['Latitude_Min']:.2f}_{row['Longitude_Min']:.2f}", axis=1)

    # Step 3: Save to a shapefile
    shapefile_name = os.path.join('shape_files',f"{element}Si_Ratio_Shpfile.shp")
    gdf.to_file(shapefile_name, driver="ESRI Shapefile")
    print(f"Shapefile saved as {shapefile_name}")

    
def create_shape_file(year,month, element,significance_level):
    # Define the path to your CSV file
    csv_file_path = f'xrf_line_catalog/{year}/{month}.csv'

    # Define the path for the output shapefile
    output_shapefile = 'Shape_Files/'+str(significance_level)+f'/{year}/{month}/{element}.shp'

    # Initialize the shapefile writer and set it to write polygons
    with shapefile.Writer(output_shapefile) as shp:
        shp.field('ID', 'N')  # Add an ID field
        shp.field('Amplitude', 'F', decimal=10)  # Add a field for Amp_data with 10 decimal places
        
        # Read rectangles from CSV and write to shapefile
        with open(csv_file_path, newline='') as csvfile:
            csv_reader = csv.reader(csvfile)
            
            # Skip header if there is one
            header = next(csv_reader)
            
            for i, row in enumerate(csv_reader):
                if i==0:
                    continue
                # Parse rectangle coordinates from CSV
                # Assuming each rectangle is defined by two points: (lat1, lon1) and (lat2, lon2)
                column_index = header.index(element)
                column_data = row[column_index]
                # print(column_index)
                                # print(column_data)
                try:
                    column_data=float(column_data)
                    Amp_index = column_index+1
                    Amp_data = float(row[Amp_index])
                    lat1,  lat2, lat3, lat4,lon1, lon2, lon3, lon4 = map(float, row[0:8])
                    
                    # Define the four corners of the rectangle
                    rectangle = [
                        [lon1, lat1],  # Bottom-left
                        [lon2, lat2],  # Bottom-right
                        [lon3, lat3],  # Top-right
                        [lon4, lat4],  # Top-left
                        [lon1, lat1]   # Close the polygon by returning to the first point
                    ]
                    
                    # Write the rectangle polygon to the shapefile
                    shp.poly([rectangle])
                    shp.record(i,Amp_data,Amp_data)  # Use the row index as ID
                except ValueError:
                    pass

    print("Shapefile created successfully.")
def combine_shapefile(category,element,significance_level):
    files=[]
    for CLASS in category:
        element_shp = CLASS+'/Shape_Files/'+str(significance_level)+'/xrf_'+CLASS+'_catalogue_'+element+'.shp'
        files.append(element_shp)

    combined_shp = 'Combined_Shape_Files/'+str(significance_level)+'/_combined_'+element+'.shp'
    with shapefile.Writer(combined_shp) as writer:
        
        # Initialize variable to check if field names are already added
        fields_added = False
        
        # Iterate over each shapefile
        for shp in files:
            # Open the shapefile
            with shapefile.Reader(shp) as reader:
                
                # Copy fields only once (first shapefile fields)
                if not fields_added:
                    writer.fields = reader.fields[1:]  # Skip the deletion flag field
                    fields_added = True
                
                # Append each shape and its corresponding record
                for shape_rec in reader.iterShapeRecords():
                    writer.shape(shape_rec.shape)  # Add the shape
                    writer.record(*shape_rec.record)  # Add the record (attributes)

    # Combined shapefile will be saved as 'combined_shapefile.shp'
    print("Shapefiles have been successfully combined.")
