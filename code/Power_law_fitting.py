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
from pyspeckit import Spectrum
import warnings
from scipy.optimize import curve_fit

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

energies = [
 1.25 ,
 1.49 ,
 1.74 ,
 6.40 ,
 3.69 ,
 4.51 ,
 5.41 ,
 0.525,
]


# Returns the sigma parameter of the initial guesses while fitting the gaussian
def find_sigma(filtered_channels, filtered_counts, el_amp, el_channel):
    half_amp = el_amp/2
    index = int((el_channel)//0.0135)
    sg = 0.125
    left = np.where(filtered_counts[0:index]<=half_amp)[0]
    right = np.where(filtered_counts[index:] <= half_amp)[0]

    if len(left) > 0 and len(right) > 0:
        left_edge = filtered_channels[left[-1]]
        right_edge = filtered_channels[right[0]]
        fwhm = abs(right_edge - left_edge)
        sigma =  fwhm/2.355
        sg=sigma
    return sg


# Gaussian Fitting, returns the fitted parameters and uncertainities corresponding to the final values
def Gauss_fit(filtered_channels,filtered_counts,guesses):
    warnings.filterwarnings("ignore")

    sp = Spectrum(data=filtered_counts, xarr=filtered_channels)

    if len(guesses):
        sp.specfit.multifit(guesses=guesses,fittype='gaussian')
        sp.plotter()
        fitted_data = sp.specfit.model
        chi_square = sp.specfit.chi2
        reduced_chi_square = chi_square / sp.specfit.dof
        fitted_params = [param.value for param in sp.specfit.parinfo]
        param_errors = [param.error for param in sp.specfit.parinfo]
        sp.specfit.plot_fit(filtered_channels)
        return fitted_params, reduced_chi_square
    else:
        return None
    
# generalised power law function
def powerlaw(x,a,c,m):
    return a*((x-c)**m)

# 
def fitted_catalogue(abs_path, year, month, output_file):
    
    date_files = os.listdir(f'{year}/{month}')
    date_files.sort()
    for i in range(len(date_files)):
        date_file =  date_files[i]
        print(date_file)
        try:

            date_filenames = pd.read_csv(f'{year}/'+ month +'/' + date_file, header = None)
            try:
                prominence_file = pd.read_csv(f'xrf_line_catalog/{year}/{month}.csv')
            except FileNotFoundError:
                print(f"File not found: xrf_line_catalog/{year}/{month}.csv")
                continue  # Skip this file and move to the next one

        except ValueError:
            continue
    
        for orbit_no in range(len(date_filenames)):

            flare_files = date_filenames[0][orbit_no].replace('\n', '').replace('[', '').replace(']', '').replace('\'', '').split()
            back_files = date_filenames[1][orbit_no].replace('\n', '').replace('[', '').replace(']', '').replace('\'', '').split()

            if len(back_files)<50:   #Atleast 50 background files to ensure sufficient statistics
                continue

            for i in range(len(back_files)):
                date_data_path  = os.path.join(abs_path,f"ch2_cla_l1_{back_files[i][11:15]}_{back_files[i][15:17]}/cla/data/calibrated/{back_files[i][11:15]}/{back_files[i][15:17]}/{back_files[i][17:19]}/")
                if i==0:
                    back_data = [fits.open(date_data_path+ back_files[i])[1].data['counts']]*len(back_files)
                else:
                    back_data[i]=   fits.open(date_data_path+ back_files[i])[1].data['counts']
        
            for i in range(len(flare_files)):
                time = flare_files[i][11:15] + '-' + flare_files[i][15:17] + '-' + flare_files[i][17:22] + ':' + flare_files[i][22:24] + ':' + flare_files[i][24:26] + '.' + flare_files[i][26:29]
                
                # Checking whether the start time of the flare file lies in the prominence_file's 'Timestamp' column or not
                row = prominence_file[prominence_file['Timestamp'] == time]

                # Initialize variables with default values at the start of the loop
                Mg_Channel, Mg_Amp = None, None
                Al_Channel, Al_Amp = None, None
                Si_Channel, Si_Amp = None, None
                Fe_Channel, Fe_Amp = None, None
                Ca_Channel, Ca_Amp = None, None
                Ti_Channel, Ti_Amp = None, None
                Cr_Channel, Cr_Amp = None, None
                O_Channel, O_Amp = None, None

                if not row['Peak_Amplitude@Mg'].isnull().values[0]: # check for non NaN value of Peak_Amplitude@Mg for background subtraction
                    if row['Peak_Amplitude@Al'].values[0]<10:
                        continue
                    date_data_path  = os.path.join(abs_path,f"ch2_cla_l1_{flare_files[i][11:15]}_{flare_files[i][15:17]}/cla/data/calibrated/{flare_files[i][11:15]}/{flare_files[i][15:17]}/{flare_files[i][17:19]}/")
                    data =  fits.open(date_data_path + flare_files[i])
                    header = data[1].header
                    channel = 13.5/1000*data[1].data['channel']
                    counts  = data[1].data['counts']
                    bg = np.sum(back_data, axis=0)/len(back_files)
            
                    selected_channels = channel[(0.5 <= channel) & (channel <= 7)]
                    selected_counts = counts[(0.5 <= channel) & (channel <= 7)]
                    bg = bg[(0.5 <= channel) & (channel <= 7)]
            
                    all_channels = np.array(selected_channels)
                    all_counts = np.array(selected_counts)

                    all_counts -= bg # subtracting the total galacric cosmic ray background from the CLASS flare file spectrum
                    
                    filtered_channels = all_channels
                    filtered_counts = all_counts
                    filtered_counts[filtered_counts<0.1]=  0.1

                    # segregate the spectrum from the peaks for accurate power law fitting
                    all_counts = np.append(all_counts[(all_channels<1)], all_counts[(all_channels>2)])
                    all_channels  = np.append(all_channels[(all_channels<1)], all_channels[(all_channels>2)])

                    Mg_Amp = row['Peak_Amplitude@Mg'].values[0]
                    Mg_Channel = row['Mg'].values[0]
                    Mg_sigma = min(find_sigma(filtered_channels, filtered_counts, Mg_Amp, Mg_Channel)/4, 0.125)
                    initial_guesses= [Mg_Amp, element_channels['Mg'], Mg_sigma]
    
                    element_exists = [False]*8
                    element_exists[0] = True

                    # limiting the fitting range to enhance accuracy and enhanced computation
                    limits = [(0,0), (element_channels['Mg'], element_channels['Mg']+0.01), (0,0.125)]
                    limited = [(True,False), (True,True), (True,True)]

                    # checking whether the element counts are not null and accordingly appending it's initial guesses for gaussian fitting later
                    if not row['Peak_Amplitude@Al'].isnull().values[0]:
                        Al_Amp = row['Peak_Amplitude@Al'].values[0]
                        Al_Channel = row['Al'].values[0]
                        Al_sigma = min(find_sigma(filtered_channels, filtered_counts, Al_Amp, Al_Channel)/4, 0.125)
                        initial_guesses_Al = [Al_Amp, element_channels['Al'], Al_sigma]
                        initial_guesses = np.concatenate([initial_guesses, initial_guesses_Al])
                        limits = np.concatenate( [limits, [(0,0), (element_channels['Al'], element_channels['Al']+0.01), (0,0.125)]], axis=0)
                        limited =np.concatenate( [limited, [(True,False), (True,True), (True,True)]], axis=0)
                        element_exists[1] = True

                    if not row['Peak_Amplitude@Si'].isnull().values[0]:
                        Si_Amp = row['Peak_Amplitude@Si'].values[0]
                        Si_Channel = row['Si'].values[0]
                        Si_sigma = min(find_sigma(filtered_channels, filtered_counts, Si_Amp, Si_Channel)/4, 0.125)
                        initial_guesses_Si = [Si_Amp, element_channels['Si'], Si_sigma]
                        initial_guesses = np.concatenate([initial_guesses, initial_guesses_Si])
                        limits = np.concatenate( [limits, [(0,0), (element_channels['Si'], element_channels['Si']+0.01), (0,0.125)]], axis=0)
                        limited =np.concatenate( [limited, [(True,False), (True,True), (True,True)]], axis=0)
                        element_exists[2] = True


                    if not row['Peak_Amplitude@Fe'].isnull().values[0]:
                        Fe_Amp = row['Peak_Amplitude@Fe'].values[0]
                        Fe_Channel = row['Fe'].values[0]
                        Fe_sigma = min(find_sigma(filtered_channels, filtered_counts, Fe_Amp, Fe_Channel)/4, 0.125)
                        initial_guesses_Fe = [Fe_Amp, element_channels['Fe'], Fe_sigma]
                        initial_guesses = np.concatenate([initial_guesses, initial_guesses_Fe])
                        limits = np.concatenate( [limits, [(0,0), (element_channels['Fe'], element_channels['Fe']+0.01), (0,0.125)]], axis=0)
                        limited =np.concatenate( [limited, [(True,False), (True,True), (True,True)]], axis=0)
                        element_exists[3] = True
                        all_counts = np.append(all_counts[(all_channels<(element_channels['Fe']-0.5))], all_counts[(all_channels>(element_channels['Fe']+0.5))])
                        all_channels  = np.append(all_channels[(all_channels<(element_channels['Fe']-0.5))], all_channels[(all_channels>(element_channels['Fe']+0.5))])


                    if not row['Peak_Amplitude@Ca'].isnull().values[0]:
                        Ca_Amp = row['Peak_Amplitude@Ca'].values[0]
                        Ca_Channel = row['Ca'].values[0]
                        Ca_sigma = min(find_sigma(filtered_channels, filtered_counts, Ca_Amp, Ca_Channel)/4, 0.125)
                        initial_guesses_Ca = [Ca_Amp, element_channels['Ca'], Ca_sigma]
                        initial_guesses = np.concatenate([initial_guesses, initial_guesses_Ca])
                        limits = np.concatenate( [limits, [(0,0), (element_channels['Ca'], element_channels['Ca']+0.01), (0,0.125)]], axis=0)
                        limited =np.concatenate( [limited, [(True,False), (True,True), (True,True)]], axis=0)
                        element_exists[4] = True
                        all_counts = np.append(all_counts[(all_channels<(element_channels['Ca']-0.5))], all_counts[(all_channels>(element_channels['Ca']+0.5))])
                        all_channels  = np.append(all_channels[(all_channels<(element_channels['Ca']-0.5))], all_channels[(all_channels>(element_channels['Ca']+0.5))])


                    if not row['Peak_Amplitude@Ti'].isnull().values[0]:
                        Ti_Amp = row['Peak_Amplitude@Ti'].values[0]
                        Ti_Channel = row['Ti'].values[0]
                        Ti_sigma = min(find_sigma(filtered_channels, filtered_counts, Ti_Amp, Ti_Channel)/4, 0.125)
                        initial_guesses_Ti = [Ti_Amp, element_channels['Ti'], Ti_sigma]
                        initial_guesses = np.concatenate([initial_guesses, initial_guesses_Ti])
                        limits = np.concatenate( [limits, [(0,0), (element_channels['Ti'], element_channels['Ti']+0.01), (0,0.125)]], axis=0)
                        limited =np.concatenate( [limited, [(True,False), (True,True), (True,True)]], axis=0)
                        element_exists[5] = True
                        all_counts = np.append(all_counts[(all_channels<(element_channels['Ti']-0.5))], all_counts[(all_channels>(element_channels['Ti']+0.5))])
                        all_channels  = np.append(all_channels[(all_channels<(element_channels['Ti']-0.5))], all_channels[(all_channels>(element_channels['Ti']+0.5))])


                    if not row['Peak_Amplitude@Cr'].isnull().values[0]:
                        Cr_Amp = row['Peak_Amplitude@Cr'].values[0]
                        Cr_Channel = row['Cr'].values[0]
                        Cr_sigma = min(find_sigma(filtered_channels, filtered_counts, Cr_Amp, Cr_Channel)/4, 0.125)
                        initial_guesses_Cr = [Cr_Amp, element_channels['Cr'], Cr_sigma]
                        initial_guesses = np.concatenate([initial_guesses, initial_guesses_Cr])
                        limits = np.concatenate( [limits, [(0,0), (element_channels['Cr'], element_channels['Cr']+0.01), (0,0.125)]], axis=0)
                        limited =np.concatenate( [limited, [(True,False), (True,True), (True,True)]], axis=0)
                        element_exists[6] = True
                        all_counts = np.append(all_counts[(all_channels<(element_channels['Cr']-0.5))], all_counts[(all_channels>(element_channels['Cr']+0.5))])
                        all_channels  = np.append(all_channels[(all_channels<(element_channels['Cr']-0.5))], all_channels[(all_channels>(element_channels['Cr']+0.5))])

                    if not row['Peak_Amplitude@O'].isnull().values[0]:
                        O_Amp = row['Peak_Amplitude@O'].values[0]
                        O_Channel = row['O'].values[0]
                        O_sigma = min(find_sigma(filtered_channels, filtered_counts, O_Amp, O_Channel)/4, 0.125)
                        initial_guesses_O = [O_Amp, element_channels['O'], O_sigma]
                        initial_guesses = np.concatenate([initial_guesses, initial_guesses_O])
                        limits = np.concatenate( [limits, [(0,0), (element_channels['O'], element_channels['O']+0.01), (0,0.125)]], axis=0)
                        limited =np.concatenate( [limited, [(True,False), (True,True), (True,True)]], axis=0)
                        element_exists[7] = True

                    initial_guesses[initial_guesses==0] = 0.01
                    # extracting out the fitted parameters and corresponding errors for power law fitting
                    [[a,c,m], err ]  = curve_fit(powerlaw, (all_channels), (all_counts),p0 = (np.max(all_counts),0.0,-1.5), maxfev = 155000, bounds=((0, -np.inf, -2.5), (4*(np.max(all_counts)), 0.5, 0)))
                    
                    # Subtracting the calculated scattered solar spectrum
                    filtered_counts-= powerlaw(filtered_channels, a,c,m)
                    filtered_counts[filtered_counts<0.1]=  0.1
                    
                    limits = [tuple(float(j) for j in i) for i in limits]
                    limited = [tuple(bool(j) for j in i) for i in limited]
                   
                    sp = Spectrum(data=filtered_counts, xarr=filtered_channels)
                    
                    # Implementing multi-gaussian fit taking into account all appended initial_guesses
                    sp.specfit.multifit(fittype = 'gaussian', guesses = initial_guesses, limits = limits, limited = limited)

                    params = [param.value for param in sp.specfit.parinfo]
                    chi_square = sp.specfit.chi2
                    reduced_chi_square = chi_square / sp.specfit.dof
                    
                    # storing the obtained final parameters element-wise (amplitude, center, sigma)
                    elements = [[None]*3]*8
                    for i in range(8):
                        if element_exists[i]:
                            index= (np.argmin(abs(initial_guesses[1::3]- energies[i])))*3
                            elements[i] =  params[index:index+3]
                                      
                    j = (np.argmin(initial_guesses[1::3]))*3
                    solar_spec_params = [a,c,m]
                    
                    if not os.path.exists(f'xrf_line_catalog_power_law/{year}'):
                        os.makedirs(f'xrf_line_catalog_power_law/{year}/')
        
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
                                        'Sigma_Mg',
                                        'Al',
                                        'Peak_Amplitude@Al',
                                        'Sigma_Al',
                                        'Si',
                                        'Peak_Amplitude@Si',
                                        'Sigma_Si',
                                        'Fe',
                                        'Peak_Amplitude@Fe',
                                        'Sigma_Fe',
                                        'Ca',
                                        'Peak_Amplitude@Ca',
                                        'Sigma_Ca',
                                        'Ti',
                                        'Peak_Amplitude@Ti',
                                        'Sigma_Ti',
                                        'Cr',
                                        'Peak_Amplitude@Cr',
                                        'Sigma_Cr',
                                        'O',
                                        'Peak_Amplitude@O',
                                        'Sigma_O',
                                        'Solar_spec_Energy',
                                        'Solar_spec_Amp',
                                        'Solar_spec_Power',
                                        'Reduced_ChiSq'
                                        ]
                                        
                        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                        if os.stat(output_file).st_size == 0:
                            writer.writeheader()
                        hdu_list_sub = data

                        writer.writerow({'Mg': elements[0][1],
                                        'Peak_Amplitude@Mg' :elements[0][0],
                                        'Sigma_Mg': elements[0][2],
                                        'Al': elements[1][1],
                                        'Peak_Amplitude@Al' :elements[1][0], 
                                        'Sigma_Al': elements[1][2],
                                        'Si': elements[2][1],
                                        'Peak_Amplitude@Si' :elements[2][0],
                                        'Sigma_Si': elements[2][2],
                                        'Fe': elements[3][1],
                                        'Peak_Amplitude@Fe' :elements[3][0],
                                        'Sigma_Fe': elements[3][2], 
                                        'Ca': elements[4][1],
                                        'Peak_Amplitude@Ca' :elements[4][0],
                                        'Sigma_Ca': elements[4][2],
                                        'Ti': elements[5][1],
                                        'Peak_Amplitude@Ti' :elements[5][0],
                                        'Sigma_Ti': elements[5][2],
                                        'Cr': elements[6][1],
                                        'Peak_Amplitude@Cr' :elements[6][0],
                                        'Sigma_Cr': elements[6][2],
                                        'O':elements[7][1],
                                        'Peak_Amplitude@O' :elements[7][0],
                                        'Sigma_O' : elements[7][2],
                                        'Solar_spec_Energy' : solar_spec_params[1],
                                        'Solar_spec_Amp': solar_spec_params[0],
                                        'Solar_spec_Power': solar_spec_params[2],
                                        'Reduced_ChiSq' : reduced_chi_square,
                                        'Latitude_1': hdu_list_sub[1].header['V0_LAT'],
                                        'Latitude_2': hdu_list_sub[1].header['V1_LAT'],
                                        'Latitude_3': hdu_list_sub[1].header['V2_LAT'],
                                        'Latitude_4': hdu_list_sub[1].header['V3_LAT'],
                                        'Longitude_1': hdu_list_sub[1].header['V0_LON'],
                                        'Longitude_2': hdu_list_sub[1].header['V1_LON'],
                                        'Longitude_3': hdu_list_sub[1].header['V2_LON'],
                                        'Longitude_4': hdu_list_sub[1].header['V3_LON'],
                                        'Timestamp': hdu_list_sub[1].header['STARTIME']})
                        

# normalizing the total background subtracted CLASS spectrum with respect to the power law background component 
def normalize_peaks(path, year,month, output_file):
    
    # load the gaussian fitted catalog for normalization
    csv_path= os.path.join(path , f'{year}/{month}.csv')
    data =  pd.read_csv(csv_path)
    amp = [0]*3
    for i in range(len(data)):
        a,c,m = data['Solar_spec_Amp'][i], data['Solar_spec_Energy'][i], data['Solar_spec_Power'][i]
        for element,j in zip(list(element_channels.keys())[0:3], range(3)):
            amp[j] = data[f'Peak_Amplitude@{element}'][i]
            channel = data[element][i]
            amp[j]/= powerlaw(channel,a,c,m) # dividing the power law component with each element's peak amplitude for attaining normalization

        if not os.path.exists(f'xrf_line_catalog_power_law_normalized/{year}'):
                os.makedirs(f'xrf_line_catalog_power_law_normalized/{year}/')
        
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
                                        'Sigma_Mg',
                                        'Al',
                                        'Peak_Amplitude@Al',
                                        'Sigma_Al',
                                        'Si',
                                        'Peak_Amplitude@Si',
                                        'Sigma_Si',
                                        'Fe',
                                        'Peak_Amplitude@Fe',
                                        'Sigma_Fe',
                                        'Ca',
                                        'Peak_Amplitude@Ca',
                                        'Sigma_Ca',
                                        'Ti',
                                        'Peak_Amplitude@Ti',
                                        'Sigma_Ti',
                                        'Cr',
                                        'Peak_Amplitude@Cr',
                                        'Sigma_Cr',
                                        'O',
                                        'Peak_Amplitude@O',
                                        'Sigma_O',
                                        'Solar_spec_Energy',
                                        'Solar_spec_Amp',
                                        'Solar_spec_Power',
                                        'Reduced_ChiSq'
                                        ]
                                        
                        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                        if os.stat(output_file).st_size == 0:
                            writer.writeheader()

                        writer.writerow({'Mg': data['Mg'][i],
                                        'Peak_Amplitude@Mg' :amp[0],
                                        'Sigma_Mg': data['Sigma_Mg'][i],
                                        'Al': data['Al'][i],
                                        'Peak_Amplitude@Al' :amp[1], 
                                        'Sigma_Al': data['Sigma_Al'][i],
                                        'Si': data['Si'][i],
                                        'Peak_Amplitude@Si' :amp[2],
                                        'Sigma_Si': data['Sigma_Si'][i],
                                        'Fe': data['Fe'][i],
                                        'Peak_Amplitude@Fe' :data['Peak_Amplitude@Fe'][i],
                                        'Sigma_Fe': data['Sigma_Fe'][i], 
                                        'Ca': data['Ca'][i],
                                        'Peak_Amplitude@Ca' : data['Peak_Amplitude@Ca'][i],
                                        'Sigma_Ca': data['Sigma_Ca'][i],
                                        'Ti': data['Ti'][i],
                                        'Peak_Amplitude@Ti' :data['Peak_Amplitude@Ti'][i],
                                        'Sigma_Ti': data['Sigma_Ti'][i],
                                        'Cr': data['Cr'][i],
                                        'Peak_Amplitude@Cr' :data['Peak_Amplitude@Cr'][i],
                                        'Sigma_Cr': data['Sigma_Cr'][i],
                                        'O':data['O'][i],
                                        'Peak_Amplitude@O' : data['Peak_Amplitude@O'][i],
                                        'Sigma_O' : data['Sigma_O'][i],
                                        'Solar_spec_Energy' : c,
                                        'Solar_spec_Amp': a,
                                        'Solar_spec_Power': m,
                                        'Reduced_ChiSq' : data['Reduced_ChiSq'][i],
                                        'Latitude_1': data['Latitude_1'][i],
                                        'Latitude_2': data['Latitude_2'][i],
                                        'Latitude_3': data['Latitude_3'][i],
                                        'Latitude_4': data['Latitude_4'][i],
                                        'Longitude_1': data['Longitude_4'][i],
                                        'Longitude_2': data['Longitude_2'][i],
                                        'Longitude_3': data['Longitude_3'][i],
                                        'Longitude_4': data['Longitude_4'][i],
                                        'Timestamp': data['Timestamp'][i]})
