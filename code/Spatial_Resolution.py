import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle
import matplotlib.image as mpimg
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from PIL import Image
import geopandas as gpd

# Moon parameters
moon_radius_km = 1737.4

# Define the directories
def avg_km_per_degree_lon(lat1, lat2):
    avg_cos = np.mean([np.cos(np.radians(lat)) for lat in np.linspace(lat1, lat2, 100)])
    return (2 * np.pi * moon_radius_km / 360) * avg_cos
def Ratio_Mapping_with_contours(year,element,desired_grid_size_km):
    increment=desired_grid_size_km/30.325
    #1 degree=30.325km
    #1km=1/30.325km
    lat_bin_size=increment##### check for what should be the bin size
    lat_min, lat_max =-90, 90
    lon_min, lon_max = -180,180
    lat_bins = np.arange(lat_min, lat_max + lat_bin_size, lat_bin_size)
    folders = year
    # Initialize an empty list to store DataFrames
    data_frames = []

    # Loop through each folder and collect CSV files
    for folder in folders:
        # Get the full path to the folder (assuming it's in the current working directory)
        folder_path = os.path.join(os.getcwd(),"xrf_line_catalog_power_law_normalized",folder)
        print(folder_path)
        # Check if the folder exists
        if os.path.exists(folder_path):
            # List all files in the folder
            for file in os.listdir(folder_path):
                # Check if the file is a CSV file
                if file.endswith(".csv"):
                    # Read the CSV file and append the DataFrame to the list
                    file_path = os.path.join(folder_path, file)  # Full path to the file
                    data_frames.append(pd.read_csv(file_path))
        else:
            print(f"Data not found for {folder}. Continuing with the available data...")
                    

    # Concatenate all DataFrames into a single DataFrame
    data = pd.concat(data_frames, ignore_index=True)

    # Display the resulting DataFrame (optional)
    data = data[abs(data['Longitude_1'] - data['Longitude_3']) < 180]
    data = data.dropna(subset=[f"Peak_Amplitude@{element}", "Peak_Amplitude@Si"])
    data = data[(data['Peak_Amplitude@Si'] >10) ]
    data[f"{element}_Si_Ratio"] = data[f"Peak_Amplitude@{element}"] / data["Peak_Amplitude@Si"]

    # Assuming your data is already loaded into a DataFrame named 'data'
    Q1 = data[f'{element}_Si_Ratio'].quantile(0.25)
    Q3 = data[f'{element}_Si_Ratio'].quantile(0.75)
    IQR = Q3 - Q1

    # Define bounds for outliers
    lower_bound = Q1 - 1.5 * IQR
    upper_bound = Q3 + 1.5 * IQR

    # Filter the DataFrame to exclude outliers
    filtered_data_iqr = data[(data[f'{element}_Si_Ratio'] >= lower_bound) & (data[f'{element}_Si_Ratio'] <= upper_bound)]

    # Now you can create a box plot with the filtered data
    data=filtered_data_iqr

    # Calculate the {element}/Si ratio

    lat_min, lat_max =-90, 90
    lon_min, lon_max = -180,180
    # Columns for latitudes and longitudes
    lat_cols = ["Latitude_1", "Latitude_2", "Latitude_3", "Latitude_4"]
    lon_cols = ["Longitude_1", "Longitude_2", "Longitude_3", "Longitude_4"]

    # Flatten lat/lon values for bounding box calculations
    latitudes = data[lat_cols].values.flatten()
    longitudes = data[lon_cols].values.flatten()

    # Create the grid
    grid_lines = []
    for i in range(len(lat_bins) - 1):
        lat1, lat2 = lat_bins[i], lat_bins[i + 1]
        avg_km_lon = avg_km_per_degree_lon(lat1, lat2)
        lon_step_size = desired_grid_size_km / avg_km_lon
        lon_bins = np.arange(lon_min, lon_max + lon_step_size, lon_step_size)
        grid_lines.append((lat1, lat2, lon_bins))
    # Initialize sum and count grids
    grid_sum = np.zeros((len(lat_bins) - 1, max(len(lon_bins) for _, _, lon_bins in grid_lines)))
    grid_count = np.zeros_like(grid_sum)
    # Assign data to the grid
    for index, row in data.iterrows():
        # Extract latitude and longitude vertices
        lat_values = row[lat_cols].values
        lon_values = row[lon_cols].values
        el_si_ratio = row[f"{element}_Si_Ratio"]

        # Calculate the min and max latitude and longitude for the area covered by the vertices
        lat_min = lat_values.min()
        lat_max = lat_values.max()
        lon_min = lon_values.min()
        lon_max = lon_values.max()

        if(abs(lon_max-lon_min)>180):
            continue

        # Iterate through each latitude grid row
        for i, (lat1, lat2, lon_bins) in enumerate(grid_lines):
            # Check if the grid cell intersects with the area defined by the min and max latitudes
            if lat1 < lat_max and lat2 > lat_min:  # Latitude overlap
                # Determine the longitude indices that overlap with the area defined by lon_min and lon_max
                lon_start_idx = np.clip(np.digitize(lon_min, lon_bins) - 1, 0, len(lon_bins) - 2)
                lon_end_idx = np.clip(np.digitize(lon_max, lon_bins) - 1, 0, len(lon_bins) - 2)

                # Iterate over the longitudinal indices to assign values to multiple grid cells
                for j in range(lon_start_idx, lon_end_idx + 1):
                    if 0 <= j < len(lon_bins) - 1:
                        # Assign the element/Si ratio to the grid cells within the latitudinal and longitudinal bounds
                        grid_sum[i, j] += el_si_ratio
                        grid_count[i, j] += 1
    # Calculate average ratio
    grid_average = np.divide(grid_sum, grid_count, out=np.zeros_like(grid_sum), where=grid_count != 0)
    grid_data = []

    # Iterate over the grid_average matrix
    for i in range(len(grid_lines)):
        lat1, lat2, lon_bins = grid_lines[i]
        for j in range(len(lon_bins) - 1):
            if grid_count[i, j] > 0:  # Only include cells with valid data
                lon1 = lon_bins[j]
                lon2 = lon_bins[j + 1]
                avg_value = grid_average[i, j]
                
                # Append the cell coordinates and average value to the list
                grid_data.append({
                    "Latitude_Min": lat1,
                    "Latitude_Max": lat2,
                    "Longitude_Min": lon1,
                    "Longitude_Max": lon2,
                    f"Average_{element}_Si_Ratio": avg_value
                })

    # Convert the list to a DataFrame
    grid_df = pd.DataFrame(grid_data)
    if not os.path.exists('ratios/'):
            os.makedirs('ratios/')
    # Save the DataFrame to a CSV file
    output_file = os.path.join(os.getcwd(),"ratios",f"grid_average_{element}_values.csv")
    print(output_file)
    grid_df.to_csv(output_file, index=False)

    background_img_path = "basemap2.jpg"  # Replace with actual image path
    background_img = Image.open(background_img_path).convert("RGB")  # Ensure the image is in RGB format
    background_img = np.array(background_img)# Replace with actual image path
    image_extent = [-180, 180, -90, 90]  # Adjust based on the background's geographic extent

    # Load the shapefile
    shapefile_path = os.path.join(os.getcwd(),"LROC_GLOBAL_MARE_180","LROC_GLOBAL_MARE_180.SHP")  # Replace with actual path
    print(shapefile_path)
    gdf = gpd.read_file(shapefile_path)

    # Ensure CRS is in EPSG:4326 (latitude/longitude)
    # if gdf.crs != "EPSG:4326":
    #     gdf = gdf.to_crs(epsg=4326)

    # Calculate polygon areas (assumes a projected CRS; for EPSG:4326, this is approximate)
    gdf['area'] = gdf.geometry.area

    # Filter polygons based on an area threshold
    area_threshold = 200 # Adjust as needed (e.g., 0 to keep all polygons)
    large_polygons = gdf[gdf['area'] > area_threshold]

    # Create the grid data (replace with your actual grid data)
    # Example grid data (replace with your real data)
    # grid_lines = [(-90 + i * 15, -90 + (i + 1) * 15, np.linspace(-180, 180, 25)) for i in range(12)]
    # grid_count = np.random.randint(0, 10, size=(12, 24))  # Example grid counts
    # grid_average = np.random.rand(12, 24)  # Example Al/Si ratios

    # Define the colormap and normalization for grid values
    norm = Normalize(vmin=lower_bound, vmax=upper_bound)  # Adjust `vmin` and `vmax` as needed
    cmap = plt.cm.turbo

    # Plotting
    fig, ax = plt.subplots(figsize=(14, 7))

    # Plot the background image
    ax.imshow(background_img, extent=image_extent, aspect='auto', zorder=1)

    # Plot hollow polygons (boundaries only) for mare regions
    large_polygons.boundary.plot(ax=ax, color='black', linewidth=1, zorder=3)

    # Add labels for polygons if 'MARE_NAME' exists
    if 'MARE_NAME' in large_polygons.columns:
        for idx, row in large_polygons.iterrows():
            centroid = row['geometry'].centroid
            ax.text(centroid.x, centroid.y, row['MARE_NAME'], fontsize=8, ha='center', color='white', zorder=4)

    # Plot the grid with grid_average values
    for i, (lat1, lat2, lon_bins) in enumerate(grid_lines):
        for j in range(len(lon_bins) - 1):
            if grid_count[i, j] > 0:  # Only plot cells with data
                rect = Rectangle(
                    (lon_bins[j], lat1),
                    lon_bins[j + 1] - lon_bins[j],
                    lat2 - lat1,
                    color=cmap(norm(grid_average[i, j])),
                    alpha=1,
                    zorder=2,  # Ensure grid cells are below the boundaries
                )
                ax.add_patch(rect)

    # Configure plot limits and labels
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xlabel("Longitude (°)", fontsize=12)
    ax.set_ylabel("Latitude (°)", fontsize=12)
    ax.set_title(f"{element}/Si Ratio Distribution on the Lunar Surface with Mare Boundaries", fontsize=15)

    # Add colorbar for the grid
    cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, label=f"{element}/Si Ratio", orientation="vertical")
    cbar.ax.tick_params(labelsize=10)

    # Add grid lines for better readability
    ax.grid(True, linestyle="--", alpha=0.3, zorder=0)

    plt.tight_layout()
    plt.show()

def Ratio_Mapping_without_contours(year,element,desired_grid_size_km):
    increment=desired_grid_size_km/30.325
    #1 degree=30.325km
    #1km=1/30.325km
    lat_bin_size=increment##### check for what should be the bin size
    lat_min, lat_max =-90, 90
    lon_min, lon_max = -180,180
    lat_bins = np.arange(lat_min, lat_max + lat_bin_size, lat_bin_size)
    folders = year
    # Initialize an empty list to store DataFrames
    data_frames = []

    # Loop through each folder and collect CSV files
    for folder in folders:
        # Get the full path to the folder (assuming it's in the current working directory)
        folder_path = os.path.join(os.getcwd(),"xrf_line_catalog_power_law_normalized",folder)
        print(folder_path)
        # Check if the folder exists
        if os.path.exists(folder_path):
            # List all files in the folder
            for file in os.listdir(folder_path):
                # Check if the file is a CSV file
                if file.endswith(".csv"):
                    # Read the CSV file and append the DataFrame to the list
                    file_path = os.path.join(folder_path, file)  # Full path to the file
                    data_frames.append(pd.read_csv(file_path))
        else:
            print(f"Data not found for {folder}. Continuing with the available data...")
                    

    # Concatenate all DataFrames into a single DataFrame
    data = pd.concat(data_frames, ignore_index=True)

    # Display the resulting DataFrame (optional)
    data = data[abs(data['Longitude_1'] - data['Longitude_3']) < 180]
    data = data.dropna(subset=[f"Peak_Amplitude@{element}", "Peak_Amplitude@Si"])

    data = data[(data['Peak_Amplitude@Si'] >10) ]
    data[f"{element}_Si_Ratio"] = data[f"Peak_Amplitude@{element}"] / data["Peak_Amplitude@Si"]

    # Assuming your data is already loaded into a DataFrame named 'data'
    Q1 = data[f'{element}_Si_Ratio'].quantile(0.25)
    Q3 = data[f'{element}_Si_Ratio'].quantile(0.75)
    IQR = Q3 - Q1

    # Define bounds for outliers
    lower_bound = Q1 - 1.5 * IQR
    upper_bound = Q3 + 1.5 * IQR

    # Filter the DataFrame to exclude outliers
    filtered_data_iqr = data[(data[f'{element}_Si_Ratio'] >= lower_bound) & (data[f'{element}_Si_Ratio'] <= upper_bound)]

    data=filtered_data_iqr

    # Calculate the {element}/Si ratio
    data[f"{element}_Si_Ratio"] = data[f"Peak_Amplitude@{element}"] / data["Peak_Amplitude@Si"]

    lat_min, lat_max =-90, 90
    lon_min, lon_max = -180,180
    # Columns for latitudes and longitudes
    lat_cols = ["Latitude_1", "Latitude_2", "Latitude_3", "Latitude_4"]
    lon_cols = ["Longitude_1", "Longitude_2", "Longitude_3", "Longitude_4"]

    # Flatten lat/lon values for bounding box calculations
    latitudes = data[lat_cols].values.flatten()
    longitudes = data[lon_cols].values.flatten()

    # Create the grid
    grid_lines = []
    for i in range(len(lat_bins) - 1):
        lat1, lat2 = lat_bins[i], lat_bins[i + 1]
        avg_km_lon = avg_km_per_degree_lon(lat1, lat2)
        lon_step_size = desired_grid_size_km / avg_km_lon
        lon_bins = np.arange(lon_min, lon_max + lon_step_size, lon_step_size)
        grid_lines.append((lat1, lat2, lon_bins))
    # Initialize sum and count grids
    grid_sum = np.zeros((len(lat_bins) - 1, max(len(lon_bins) for _, _, lon_bins in grid_lines)))
    grid_count = np.zeros_like(grid_sum)
    # Assign data to the grid
    for index, row in data.iterrows():
        # Extract latitude and longitude vertices
        lat_values = row[lat_cols].values
        lon_values = row[lon_cols].values
        el_si_ratio = row[f"{element}_Si_Ratio"]

        # Calculate the min and max latitude and longitude for the area covered by the vertices
        lat_min = lat_values.min()
        lat_max = lat_values.max()
        lon_min = lon_values.min()
        lon_max = lon_values.max()

        if(abs(lon_max-lon_min)>180):
            continue

        # Iterate through each latitude grid row
        for i, (lat1, lat2, lon_bins) in enumerate(grid_lines):
            # Check if the grid cell intersects with the area defined by the min and max latitudes
            if lat1 < lat_max and lat2 > lat_min:  # Latitude overlap
                # Determine the longitude indices that overlap with the area defined by lon_min and lon_max
                lon_start_idx = np.clip(np.digitize(lon_min, lon_bins) - 1, 0, len(lon_bins) - 2)
                lon_end_idx = np.clip(np.digitize(lon_max, lon_bins) - 1, 0, len(lon_bins) - 2)

                # Iterate over the longitudinal indices to assign values to multiple grid cells
                for j in range(lon_start_idx, lon_end_idx + 1):
                    if 0 <= j < len(lon_bins) - 1:
                        # Assign the element/Si ratio to the grid cells within the latitudinal and longitudinal bounds
                        grid_sum[i, j] += el_si_ratio
                        grid_count[i, j] += 1
    # Calculate average ratio
    grid_average = np.divide(grid_sum, grid_count, out=np.zeros_like(grid_sum), where=grid_count != 0)

    grid_data = []

    # Iterate over the grid_average matrix
    for i in range(len(grid_lines)):
        lat1, lat2, lon_bins = grid_lines[i]
        for j in range(len(lon_bins) - 1):
            if grid_count[i, j] > 0:  # Only include cells with valid data
                lon1 = lon_bins[j]
                lon2 = lon_bins[j + 1]
                avg_value = grid_average[i, j]
                
                # Append the cell coordinates and average value to the list
                grid_data.append({
                    "Latitude_Min": lat1,
                    "Latitude_Max": lat2,
                    "Longitude_Min": lon1,
                    "Longitude_Max": lon2,
                    f"Average_{element}_Si_Ratio": avg_value
                })

    # Convert the list to a DataFrame
    grid_df = pd.DataFrame(grid_data)

    # Save the DataFrame to a CSV file
    output_file = os.path.join("ratios",f"grid_average_{element}_values.csv")
    grid_df.to_csv(output_file, index=False)


    image_path = "basemap2.jpg"
    moon_map = Image.open(image_path)
    # moon_map = moon_map.resize((moon_map.width // 2, moon_map.height // 2))  # Resize for performance
    moon_map = moon_map.convert("RGB")

    # Step 3: Convert the image to a NumPy array
    moon_map_array = np.array(moon_map)

    # Check the shape and dtype of the array
    print(f"Image shape: {moon_map_array.shape}, dtype: {moon_map_array.dtype}")
    moon_map_array = np.array(moon_map)
    fig, ax = plt.subplots(figsize=(14, 7))

    # Display the image
    ax.imshow(moon_map, extent=[-180, 180, -90, 90], aspect='auto')

    # Normalize grid values for coloring
    norm = plt.Normalize(vmin=lower_bound, vmax=upper_bound)
    cmap = plt.cm.turbo

    for i, (lat1, lat2, lon_bins) in enumerate(grid_lines):
        for j in range(len(lon_bins) - 1):
            if grid_count[i, j] > 0:  # Only plot non-empty cells
                rect = Rectangle(
                    (lon_bins[j], lat1),
                    lon_bins[j + 1] - lon_bins[j],
                    lat2 - lat1,
                    color=cmap(norm(grid_average[i, j])),
                    alpha=1,
                )
                ax.add_patch(rect)


    # Configure the plot
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xlabel("Longitude (°)")
    ax.set_ylabel("Latitude (°)")
    ax.set_title(f"{element}/Si Ratio Distribution on the Lunar Surface")
    plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, label="Al/Si Ratio")
    plt.grid(True, linestyle="--", alpha=0.3)

    # Show the plot
    plt.show()


