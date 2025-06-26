import numpy as np
import pandas as pd
from astropy.io import fits
import matplotlib.pyplot as plt
import datetime as dt
import os
from scipy.signal import find_peaks

def extract_orbits(data_path, back_angle,  year, month , day):
    date= dt.datetime(year, month, day)
    previous_date = date - dt.timedelta(1)
    next_date =  date+ dt.timedelta(1)

    date_dir = os.path.join(data_path,f'ch2_cla_l1_{date.year}_{date.month}/cla/data/calibrated/{date.year}/{date.month}/{date.day}')
    next_date_dir = os.path.join(data_path, f'ch2_cla_l1_{next_date.year}_{next_date.month}/cla/data/calibrated/{next_date.year}/{next_date.month}/{next_date.day}')
    
    # Adding a '0' to the start of the date/month given it is less than 10 eg:- 6 -> 06
    if date.day<10:
         date_dir = date_dir[0:-1] + f'0{date.day}'
    if next_date.day<10:
         next_date_dir = next_date_dir[0:-1] + f'0{next_date.day}'

    if date.month<10:

         date_dir  = date_dir[0:-31] + f'0{date.month}' + date_dir[-30:-4] + f'0{date.month}'+ date_dir[-3:]

    if next_date.month<10:
         next_date_dir  = next_date_dir[0:-31] + f'0{next_date.month}' + next_date_dir[-30:-4] + f'0{next_date.month}'+ next_date_dir[-3:]
    
     # loading file name of that particluar date
    date_filenames =  [i for i in os.listdir(date_dir) if i.endswith('fits')]
    date_filenames.sort()
    orbits = []
    angles=  [0]*len(date_filenames)        
    for i in range(len(date_filenames)):
        file = date_filenames[i]
        header = fits.open(date_dir + '/'+ file)[1].header
        angles[i] = header['phaseang'] # loading the solar angle for each file

    # Observing the periodic behaviour of the solar angle going from a particular maxima to a minima and back to that same maxima marking completion of one orbit
    peaks = find_peaks(angles)[0] # peak here refers to the maxima of the solar angle of that particular orbit
    if len(peaks)==0:
         return orbits
    for i in range(len(peaks)):
        if i!= (len(peaks)-1):
            angles_orbit = angles[peaks[i]:peaks[i+1]]
            date_filenames_orbit = np.array(date_filenames[peaks[i]:peaks[i+1]])
            back_files = (date_filenames_orbit[np.array(angles_orbit)>back_angle])
            sun_files = (date_filenames_orbit[np.array(angles_orbit)<90])
            orbits.append([sun_files, back_files])

          # observing the last peak of the day
        else:
            angles_orbit = angles[peaks[i]:-1]
            date_filenames_orbit = np.array(date_filenames[peaks[i]:-1])
            back_files = (date_filenames_orbit[np.array(angles_orbit)>back_angle]) # files considered for background
            sun_files = (date_filenames_orbit[np.array(angles_orbit)<90]) # files considered for daytime

          # checking whether the next day exists or not
            if os.path.isdir(next_date_dir):

                next_date_filenames =  [j for j in os.listdir(next_date_dir) if j.endswith('fits')]
                next_date_filenames.sort()
                angles=  [0]*len(next_date_filenames)        
                for k in range(len(next_date_filenames)):
                        file = next_date_filenames[k]
                        header = fits.open(next_date_dir + '/'+ file)[1].header
                        angles[k] = header['phaseang']
                
                peaks = find_peaks(angles)[0]
                if len(peaks)==0:
                     orbits.append([sun_files, back_files])
                     return orbits
                
                
                angles_orbit = angles[0:peaks[0]]
                next_date_filenames_orbit = np.array(next_date_filenames[0:peaks[0]])
                
                # appending the background and flare files to the same orbit
                back_files = (np.append(back_files, next_date_filenames_orbit[np.array(angles_orbit)>back_angle]))
                sun_files = (np.append(sun_files, next_date_filenames_orbit[np.array(angles_orbit)<90]))

               # appending the orbit to list of all orbits
                orbits.append([sun_files, back_files])

            else:
                 
                 # appending the orbit to list of all orbits if next day not found
                 orbits.append([sun_files, back_files]) 


    return orbits