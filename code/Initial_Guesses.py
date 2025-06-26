from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os
from datetime import datetime, timedelta
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.stats import tstd
import csv
from tqdm import tqdm
from scipy.signal import find_peaks
from Create_shp import *

element_channels = {
    'Mg': (1.25 ),  # 1.25 keV
    'Al': (1.49 ),  # 1.49 keV
    'Si': (1.74 ),  # 1.74 keV
    'Ca': (3.69 ),  # 3.69 keV
    'Ti': (4.51 ),  # 4.51 keV
    'Cr': (5.41 ),  # 5.41 keV
    'Fe': (6.40 ),  # 6.40 keV
    'O': (0.525 )   # 0.525 keV
}

def find_nearest_peak(final_peaks,final_channels,element_channels,element):
    target_energy=element_channels[element]
        # Find the nearest peak to the target energy
    if len(final_peaks) > 0:
        selected_channels =  abs(final_channels-target_energy)<=0.125
        final_channels = final_channels[selected_channels]
        final_peaks  = final_peaks[selected_channels]

        if final_channels.size>0:
             return final_peaks[np.argmin(abs(final_channels-target_energy)/final_peaks)], final_channels[np.argmin(abs(final_channels-target_energy)/final_peaks)]
        else:
            return None,None
    else:
        return None, None
    
    
def generate_csv(filename):
    with open(filename, mode='w', newline='') as csvfile:
        pass

def coverage_catalogue(abs_path,year, month, output_file,significance_level):    
    
    # Listing out the CLASS files for a given year and month
    date_files = os.listdir(f'{year}/{month}')

    # Sorting the date_files for analysis
    date_files.sort()
    for i in range(len(date_files)):
        date_file =  date_files[i] #Using 1 file at a time
        print(date_file)
        try:

            date_filenames = pd.read_csv(f'{year}/'+ month +'/' + date_file, header = None) #Loading the dataframe

        except ValueError:
            continue

        # 
        for orbit_no in range(len(date_filenames)):

            #Segregating the flare and the background files orbit-wise
            flare_files = date_filenames[0][orbit_no].replace('\n', '').replace('[', '').replace(']', '').replace('\'', '').split()
            back_files = date_filenames[1][orbit_no].replace('\n', '').replace('[', '').replace(']', '').replace('\'', '').split()

            if len(back_files)<50:   #Atleast 50 background files to ensure sufficient statistics
                continue

            # Effective cosmic ray background for the particular orbit
            for i in range(len(back_files)):
                date_data_path  = os.path.join(abs_path,f"ch2_cla_l1_{back_files[i][11:15]}_{back_files[i][15:17]}/cla/data/calibrated/{back_files[i][11:15]}/{back_files[i][15:17]}/{back_files[i][17:19]}/")
                if i==0:
                    back_data = [fits.open(date_data_path+ back_files[i])[1].data['counts']]*len(back_files)
                else:
                    back_data[i]=   fits.open(date_data_path+ back_files[i])[1].data['counts']
            
            # Processing the flare files based on the background continuum obtained previously.
            for i in range(len(flare_files)):
                date_data_path  = os.path.join(abs_path,f"ch2_cla_l1_{flare_files[i][11:15]}_{flare_files[i][15:17]}/cla/data/calibrated/{flare_files[i][11:15]}/{flare_files[i][15:17]}/{flare_files[i][17:19]}/")
                data =  fits.open(date_data_path+ flare_files[i])
                header = data[1].header
                channel = 13.5/1000*data[1].data['channel']
                counts  = data[1].data['counts']
                bg = significance_level*tstd(back_data)+np.sum(back_data, axis=0)/len(back_files)
        
                selected_channels = channel[(0.5 <= channel) & (channel <= 7)]
                selected_counts = counts[(0.5 <= channel) & (channel <= 7)]
                bg=  bg[(0.5 <= channel) & (channel <= 7)]
        
        
                all_channels = np.array(selected_channels)
                all_counts = np.array(selected_counts)

                # Ensuring dynamic peak amplitude handling of the CLASS spectrum
                peak_indices = find_peaks(all_counts, prominence=all_counts.mean())[0]
                filtered_channels = all_channels[peak_indices]
                filtered_counts = all_counts[peak_indices]

                Mg_Amp,Mg_Channel=find_nearest_peak(filtered_counts ,filtered_channels, element_channels,'Mg') 
                Al_Amp,Al_Channel=find_nearest_peak( filtered_counts,filtered_channels,element_channels,'Al')
                Si_Amp,Si_Channel=find_nearest_peak( filtered_counts,filtered_channels,element_channels,'Si')           
                Ca_Amp,Ca_Channel=find_nearest_peak( filtered_counts,filtered_channels,element_channels,'Ca')
                Ti_Amp,Ti_Channel=find_nearest_peak( filtered_counts,filtered_channels,element_channels,'Ti')
                Cr_Amp,Cr_Channel=find_nearest_peak( filtered_counts,filtered_channels,element_channels,'Cr')
                Fe_Amp,Fe_Channel=find_nearest_peak( filtered_counts,filtered_channels,element_channels,'Fe')
                O_Amp,O_Channel=find_nearest_peak( filtered_counts,filtered_channels,element_channels,'O')

                if not os.path.exists(f'xrf_line_catalog/{year}'):
                    os.makedirs(f'xrf_line_catalog/{year}/')
    
                with open(output_file, 'a', newline='') as csvfile:
                    fieldnames = ['Latitude_1',
                                    'Latitude_2',
                                    'Latitude_3',
                                    'Latitude_4',
                                    'Longitude_1',
                                    'Longitude_2',
                                    'Longitude_3',
                                    'Longitude_4' ,
                                    'Timestamp',
                                    'Mg',
                                    'Peak_Amplitude@Mg',
                                    'Al',
                                    'Peak_Amplitude@Al',
                                    'Si',
                                    'Peak_Amplitude@Si',
                                    'Fe',
                                    'Peak_Amplitude@Fe',
                                    'Ca',
                                    'Peak_Amplitude@Ca',
                                    'Ti',
                                    'Peak_Amplitude@Ti',
                                    'Cr',
                                    'Peak_Amplitude@Cr',
                                    'O',
                                    'Peak_Amplitude@O',
                                    ]
                                    
                    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                    if os.stat(output_file).st_size == 0:
                        writer.writeheader()
                    hdu_list_sub=data

                    # Storing Background Subtracted Peaks
                    if bg[Mg_Channel==selected_channels].size==0:
                        pass
                    elif Mg_Amp <bg[Mg_Channel==selected_channels][0]:
                        Mg_Amp, Mg_Channel = [None,None]
                    elif Mg_Amp >= bg[Mg_Channel==selected_channels][0]:
                        Mg_Amp-= (bg[Mg_Channel==selected_channels])[0]

                    if bg[Al_Channel==selected_channels].size==0:
                        pass
                    elif Al_Amp < bg[Al_Channel==selected_channels][0]:
                        Al_Amp,Al_Channel = [None,None]
                    elif Al_Amp >= bg[Al_Channel==selected_channels][0]:
                        Al_Amp-= (bg[Al_Channel==selected_channels])[0]

                    if bg[Si_Channel==selected_channels].size==0:
                        pass
                    elif Si_Amp < bg[Si_Channel==selected_channels][0]:
                        Si_Amp,Si_Channel =  [None,None]
                    elif Si_Amp >= bg[Si_Channel==selected_channels][0]:
                        Si_Amp-= bg[Si_Channel==selected_channels][0]

                    if bg[Fe_Channel==selected_channels].size==0:
                        pass
                    elif Fe_Amp < bg[Fe_Channel==selected_channels][0]:
                        Fe_Amp,Fe_Channel=  [None,None]
                    elif Fe_Amp >= bg[Fe_Channel==selected_channels][0]:
                        Fe_Amp-= bg[Fe_Channel==selected_channels][0]

                    if bg[Ca_Channel==selected_channels].size==0:
                        pass
                    elif Ca_Amp < bg[Ca_Channel==selected_channels][0]:
                        Ca_Amp,Ca_Channel =  [None,None]
                    elif Ca_Amp > bg[Ca_Channel==selected_channels][0]:
                        Ca_Amp-=bg[Ca_Channel==selected_channels][0]
                    
                    if bg[Ti_Channel==selected_channels].size==0:
                        pass
                    elif Ti_Amp < bg[Ti_Channel==selected_channels][0]:
                        Ti_Amp,Ti_Channel= [None,None]
                    elif Ti_Amp >= bg[Ti_Channel==selected_channels][0]:
                        Ti_Amp-=bg[Ti_Channel==selected_channels][0]
                    
                    if bg[O_Channel==selected_channels].size==0:
                        pass
                    elif O_Amp < bg[O_Channel==selected_channels][0]: 
                        O_Amp,O_Channel= [None,None]
                    elif O_Amp >= bg[O_Channel==selected_channels][0]:
                        O_Amp-= bg[O_Channel==selected_channels][0]
                        
                    if bg[Cr_Channel==selected_channels].size==0:
                        pass
                    elif Cr_Amp < bg[Cr_Channel==selected_channels][0]:
                        Cr_Amp,Cr_Channel= [None,None]
                    elif Cr_Amp >= bg[Cr_Channel==selected_channels][0]:
                        Cr_Amp-= bg[Cr_Channel==selected_channels][0]
                    
                    if any(v is None for v in [Mg_Amp, Al_Amp, Si_Amp]):
                        Mg_Amp, Mg_Channel = [None,None]
                        Al_Amp,Al_Channel = [None,None]
                        Si_Amp,Si_Channel =  [None,None]
                        Fe_Amp,Fe_Channel=  [None,None]
                        Ca_Amp,Ca_Channel =  [None,None]
                        Ti_Amp,Ti_Channel= [None,None]
                        O_Amp,O_Channel= [None,None]
                        Cr_Amp,Cr_Channel= [None,None]
                    else:
                        pass
                        
                    writer.writerow({'Mg': Mg_Channel,
                                    'Peak_Amplitude@Mg': Mg_Amp,
                                    'Al': Al_Channel,
                                    'Peak_Amplitude@Al': Al_Amp, 
                                    'Si': Si_Channel,
                                    'Peak_Amplitude@Si': Si_Amp,
                                    'Fe': Fe_Channel,
                                    'Peak_Amplitude@Fe': Fe_Amp, 
                                    'Ca': Ca_Channel,
                                    'Peak_Amplitude@Ca': Ca_Amp,
                                    'Ti': Ti_Channel,
                                    'Peak_Amplitude@Ti': Ti_Amp,
                                    'Cr': Cr_Channel,
                                    'Peak_Amplitude@Cr': Cr_Amp,
                                    'O': O_Channel,
                                    'Peak_Amplitude@O': O_Amp,
                                    'Latitude_1': hdu_list_sub[1].header['V0_LAT'],
                                    'Latitude_2': hdu_list_sub[1].header['V1_LAT'],
                                    'Latitude_3': hdu_list_sub[1].header['V2_LAT'],
                                    'Latitude_4': hdu_list_sub[1].header['V3_LAT'],
                                    'Longitude_1': hdu_list_sub[1].header['V0_LON'],
                                    'Longitude_2': hdu_list_sub[1].header['V1_LON'],
                                    'Longitude_3': hdu_list_sub[1].header['V2_LON'],
                                    'Longitude_4': hdu_list_sub[1].header['V3_LON'],
                                    'Timestamp': hdu_list_sub[1].header['STARTIME']})