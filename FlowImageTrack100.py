#!/usr/bin/env python
import numpy as np
import math as mm
from scipy import stats

def main():
   """Sort flow imaging files for sedimentation or 100% counting analysis.

    PROGRAM: FlowImageTrack100

    VERSION: 1.1
    DATE:  08 JAN 2018
    PURPOSE:  Sorting of flow imaging exported data files for analyzing apparent
    particle densities by sedimentation properties and for performing 100%
    particle counts
    INPUT & OUTPUT FILES:  Both are .csv format
    PYTHON VERSION:  Tested on Python 3.6
    INPUT FILE: specified number of header rows; succeeding rows have:
      particle image id, area (um^2), x & y corner location (pixels),diameter
      (um), elapsed time (seconds), image height & width (pixels) 
    OUTPUT FILE: 
      Input file name echoed
      Allowable change in x & y between images (pixels), allowable fractional
      change in dia., time to look ahead for next matching image (seconds)
      Rows reporting: Part. #, ID 1st image, # images in track, ave. time
      (seconds), ave. x (pixels), ave. diameter (um), slope of linear fit
      (um per second), rms dev. from fit
    Optional output section containing: 
      ID, diameter (um), area-based diameter (ABD) (um), area (um^2), center
      x (pixels), time (seconds), y (pixels), x_corner (pixels), y_corner
      (pixels), track number

      For 100% counting, the number of particles found = the total number of
      tracks found that have at least the specified n_min number of tracks.
      Since the track number will increment even if there is a single particle
      in the track, the total particle number is not equal to the highest track
      number.

    2/14/2017  Code translated from Fortran FlowCamSed.for
    AUTHOR:  Dean Ripple, NIST

    Names are used consistently with certain Types; 
    square brackets indicate Names used as Lists:
    String: filein, fileout, ctrack
    Integer: mask[], id[], IDcounter[], IDparticle[], ID0[]
    Float: dia[], area[], x[], y[], xt[], yt[], ave_dia, Slope, time[],
      rms, ave_time, ave_x, dia_abd[], xcorner[], ycorner[]
    Boolean: last_image, done, done2
    """
   
pi = 3.14159
n_min = 2    #Minimum number of images to count a particle
n_delta_min = 10  #Minimum change in y pixels to count a particle
n_y_min = 50  #Minimum y pixel to count an image
n_y_max = 900 #Maximum y pixel to count an image

image_id = []
area = []
xcorner = []
ycorner = []
x = []
y = []
dia = []
dia_abd = []
time = []
mask = []
xt = []
yt = []
id0 = []
filein = input('Enter file name for data: ')
num_header = int(input("Enter number of header lines: "))
fileino = open(filein, "r")
fileout = input('Enter output file name: ')
# Read header lines from input data file:
for i in range(num_header):
    fileino.readline()
# Now read in numerical data
for line in fileino:
    lineList = line.split(",")
# Calculate middle y pixel value of image, & proceed if it's in range
    ywidth = float(lineList[6])
    ytemp = float(lineList[3]) + ywidth/2.
    if ytemp >= n_y_min and ytemp <= n_y_max:
        image_id.append(int(lineList[0]))
        n_count = len(image_id)
        area.append(float(lineList[1]))
        xcorner.append(float(lineList[2]))
        ycorner.append(float(lineList[3]))
        dia.append(float(lineList[4]))
        dia_abd.append(2.*mm.sqrt(area[n_count - 1]/pi))
        time.append(float(lineList[5]))
        xwidth = float(lineList[7])
        x.append(xcorner[n_count - 1] + xwidth/2.)
        y.append(ytemp)
        mask.append(int(0))
#        print('{0:4d}, {1:6.2f}, {2:6.1f}, {3:6.1f}'.format(
#              image_id[n_count - 1], area[n_count -1], x[n_count - 1],
#              ytemp))     
fileino.close()
# blank element on time[] & mask[] to avoid an end-of-list crash below
time.append(0.0)
mask.append(0)
time_seek = float(input('Number of seconds to look ahead for next image: '))
n_delta_x = int(input('Threshold in pixels for x variations: '))
n_delta_y = int(input('Threshold in pixels for y variations: '))
# diameters allowed in range 1 - delta_d <(present diameter)/(prior diameter)
#  1 + delta_d
delta_d = float(input('Threshold for diameter variation (e.g., 0.2): '))
particle_num = 1
index = 0
done = False
while not done:
# Identify next unsorted particle
    id0.append(image_id[index])    
    dia0 = dia[index]
    x0 = x[index] 
    y0 = y[index]
    mask[index] = particle_num
##    print('Particle # & 1st image ID: {0:4d}, {1:4d}'.
##      format(particle_num, id0[particle_num-1]))
# Tag all particles within next time_seek that match x, y, & diameter of previous image:
    last_image =  False
    last_i = index
    while not last_image:
        last_image = True
        i = last_i
        while (time[i] <= (time[last_i] + time_seek)) and (i <= n_count - 1):
            dia_ration = dia[i]/dia0
            if (mask[i] == 0 and abs(x[i] - x0) < n_delta_x and 
                dia_ration < (1 + delta_d) and dia_ration > (1 - delta_d) and
                (y[i] - y0) < n_delta_y and (y[i] - y0) > -n_delta_x):
                mask[i] = particle_num
                x0 = x[i]
                y0 = y[i]
                dia0 = dia[i]
                last_i = i
                last_image = False
                break
            else:
                i += 1 
# Find next available index=image count
    done2 = False
    while mask[index] != 0 and not done2:            
        index += 1
        if index > n_count - 1:
            print('Ran out of points at Particle # = {0:4d}'.
                  format(particle_num))
            done2 = True
            done = True        
    if not done: 
        particle_num += 1
n_particle_count = particle_num
# Now loop through particles, and compute average diameters and average velocity
# (y vs. time)       
fileouto = open(fileout, 'a')
#print('Part. #, 1st image ID, # images, ave. time, ave. x, ave. dia., slope, rms residuals')
print('Working!')
fileouto.write('Input file name: {0}\n'.format(filein))
fileouto.write('Delta x, Delta y, Delta d, time to look ahead: \n')
fileouto.write('{0:3d}, {1:3d}, {2:7.3f}, {3:7.3f}\n'.format(
               n_delta_x, n_delta_y, delta_d, time_seek))
fileouto.write('Part. #, ID 1st image, # images, ave. time, ave. x, ' + 
               'ave. dia., slope, rms residuals\n')
n_full_tracks = 0
for i in range(n_particle_count+1):
    num_set = 0
    ave_dia_sum = 0.
    ave_x_sum = 0.
    ave_time_sum = 0.
    min_y = 2000
    max_y = 0
    xt = []
    yt = []
    for j in range(n_count):
        if mask[j] == i:
            num_set +=  1
            ave_dia_sum = ave_dia_sum + dia[j]
            ave_x_sum = ave_x_sum + x[j]
            ave_time_sum = ave_time_sum + time[j]
            xt.append(time[j])
            yt.append(y[j])
            if y[j] < min_y: 
                min_y = y[j]
            if y[j] > max_y: 
                max_y = y[j]
    n_delta_y = max_y - min_y
# Require at least n_min points over a span of n_delta_min pixels to count
    if num_set >= n_min and n_delta_y > n_delta_min:
        ave_dia = ave_dia_sum/num_set
        ave_time = ave_time_sum/num_set
        ave_x = ave_x_sum/num_set
        n_full_tracks += 1
# Least squares fit of slope, if there are enough points 
        if num_set >= 3:
            slope, intercept, r_value, p_value, rms = stats.linregress(xt,yt)
        else:
            if xt[1] != xt[0]:
               slope = (yt[1] - yt[0])/(xt[1] - xt[0])
            else:
               slope = 1e9
            rms = 1e10 
        msg = ('{0:6d}, {1:8d}, {2:8d}, {3:12.4g}, {4:7.2f}, {5:7.2f},' 
               '{6:12.3g}, {7:14.3g},'.format(i, id0[i - 1], num_set, ave_time,
               ave_x, ave_dia, slope, rms))      
#        print(msg)
        msg += '\n'
        fileouto.write(msg)
msg = 'Number of tracks with fits: ' + '{:d}'.format(n_full_tracks)
msg += '\n'
fileouto.write(msg)
print(msg)
# Echo input data, if desired                  
ctrack = input('Write data set with track ID as last column (y/n)?')
if ctrack == 'y':
    fileouto.write('Image ID, dia, ABD dia, area, time, x, y, x_corner, ' +
                   'y_corner, track #\n')
    for i in range(n_count):
        msg = ('{0:8d}, {1:7.2f}, {2:7.2f}, {3:12.4g}, {4:8.2f}, {5:8.2f},'
               ' {6:8.2f}, {7:8.2f}, {8:8.2f}, {9:5d},\n'.format(
               image_id[i], dia[i], dia_abd[i], area[i], time[i], x[i], y[i], 
               xcorner[i], ycorner[i], mask[i])) 
        fileouto.write(msg)   
fileouto.close()

if __name__ == "__main__":
    main()
