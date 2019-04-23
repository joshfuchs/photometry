# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 11:18:53 2018

@author: rrolvera
"""

import numpy as np # You always have to do this, the np is to avoid confusion
import os # You also have to do this, it stands for operating system
from astropy.io import fits # Imported so we can open the fits files
from photutils import centroid_2dg, CircularAperture, CircularAnnulus # for the centroid
import matplotlib.pyplot as plt 
from photutils import aperture_photometry
import math 
from astropy.stats import LombScargle


print("Part 1:")

# First, we want to summon the list that we want, once we find it, we want to print it
# "Read text file of images"

print("Step 1:")

def Read_List( lst ):
    # This function reads a list of images and decomposes them into a python list of image names - meza
    # This works as a function because it is doing as the command says 
    list_file = open(lst,'r') # This opens the file 
    im_list = list_file.read() #This reads the file
    list_file.close() # This closes the file
    im_list = im_list.split() # This makes the file into a list of some sort
    return im_list # This produces the list we just found 

os.chdir('/home/stdres1/AstroResearch/WD1150-153/photometry') # Changes our directory to where our "listphot" is located
images = 'listphot' # A name is given to the file we are using/want 
image_list = Read_List(images) # Read_List command reads dummy and is called dummy_lists as a whole
print(image_list) # This reads out our list

# Second, we want to define the four stars we will be using
# "Set coordinates of initial target + 3 comparisons"

print("Step 2:")

COORDINATES = 'coordinates.txt'
x1, y1, x2, y2, x3, y3, x4, y4 = np.genfromtxt(COORDINATES, dtype=int, comments='#', unpack=True)
target_star = [] # This lets us set a name to a set of coordinates
target_star.append((x1, y1)) # This adds our coordinate/star to our program
print("Target Star Coordinates:") # This is the title of information
print(target_star) # This prints our star coordinates

comparison_star_one = [] # This lets us set a name to a set of coordinates
comparison_star_one.append((x2, y2)) # This adds our coordinate/star to our program
print("Comparison Star One Coordinates:") # This is the title of information
print(comparison_star_one) # This prints our star coordinates
        
comparison_star_two = [] # This lets us set a name to a set of coordinates
comparison_star_two.append((x3, y3)) # This adds our coordinate/star to our program
print("Comparison Star Two Coordinates:") # This is the title of information
print(comparison_star_two) # This prints our star coordinates

comparison_star_three = [] # This lets us set a name to a set of coordinates
comparison_star_three.append((x4, y4)) # This adds our coordinate/star to our program
print("Comparison Star Three Coordinates:") # This is the title of information
print(comparison_star_three) # This prints our star coordinates        

# Third, we want to create some lists

print("Step 3:")

time = [] # This initilizes the array, in the future we can input information

target_flux = [] # This initilizes the array, in the future we can input information
c1_flux = []
c2_flux = []
c3_flux = []

target_sky = []
c1_sky = []
c2_sky = []
c3_sky = []

actual_target = []
actual_c1 = []
actual_c2 = []
actual_c3 = []

sky_counts_target_lst = []
sky_counts_c1_lst = []
sky_counts_c2_lst = []
sky_counts_c3_lst= []


# Now we can write a "for" loop 

print("Part 2:")

# First, we want to write a "for" loop that reads out the list as statements
# Read in Fits File (astropy.io.fits)
# How in the universe to do insert open.fits into this loop?

print("Step 1:")
pix_dis = 7
annu_radi_1 = 30
annu_radi_2 = 35
images = image_list # This defines that the word dummy is the dummy_lists
for words in images: # This "for" loop is for the statements in dummy_lists
    hdu = fits.getheader(words) # This gets the words from the header
    img_data = fits.getdata(words) # This reads it as data
    print('=====================================Newest Update=====================================\n')
    print('Image Dimensions:')
    print(img_data.shape) # This shows us the the list
    print('\n********************Centroiding********************\n')
    marker = '+'
    ms, mew = 30, 2.
    
    targetx1, targety1 = centroid_2dg(img_data[0,y1-pix_dis:y1+pix_dis,x1-pix_dis:x1+pix_dis])
    print('Centroid Target Star:')
    fig, ax = plt.subplots(1, 1)
    ax.imshow(img_data[0,y1-pix_dis:y1+pix_dis,x1-pix_dis:x1+pix_dis], origin='lower', interpolation='nearest', cmap='viridis')
    plt.plot(targetx1, targety1, color='#17becf', marker=marker, ms=ms, mew=mew)
    plt.show()
    print(targetx1, targety1)
    
    targetx2, targety2 = centroid_2dg(img_data[0,y2-pix_dis:y2+pix_dis,x2-pix_dis:x2+pix_dis])
    print('\nCentroid Comparison Star 1:')
    fig, ax = plt.subplots(1, 1)
    ax.imshow(img_data[0,y2-pix_dis:y2+pix_dis,x2-pix_dis:x2+pix_dis], origin='lower', interpolation='nearest', cmap='viridis')
    plt.plot(targetx2, targety2, color='r', marker=marker, ms=ms, mew=mew)
    plt.show()
    print(targetx2, targety2)
    
    targetx3, targety3 = centroid_2dg(img_data[0,y3-pix_dis:y3+pix_dis,x3-pix_dis:x3+pix_dis])
    print('\nCentroid Comparison Star 2:')
    fig, ax = plt.subplots(1, 1)
    ax.imshow(img_data[0,y3-pix_dis:y3+pix_dis,x3-pix_dis:x3+pix_dis], origin='lower', interpolation='nearest', cmap='viridis')
    plt.plot(targetx3, targety3, color='y', marker=marker, ms=ms, mew=mew)
    plt.show()
    print(targetx3, targety3)
    
    targetx4, targety4 = centroid_2dg(img_data[0,y4-pix_dis:y4+pix_dis,x4-pix_dis:x4+pix_dis])
    print('\nCentroid Comparison Star 3:')
    fig, ax = plt.subplots(1, 1)
    ax.imshow(img_data[0,y4-pix_dis:y4+pix_dis,x4-pix_dis:x4+pix_dis], origin='lower', interpolation='nearest', cmap='viridis')
    plt.plot(targetx4, targety4, color='k', marker=marker, ms=ms, mew=mew)
    plt.show()
    print(targetx4, targety4)
    
    positions = [(x1 + targetx1, y1 + targety1)]
    positions_c1 =  [(x2 + targetx2, y2 + targety2)]
    positions_c2 = [(x3 + targetx3, y3 + targety3)]
    positions_c3 = [(x4 + targetx4, y4 + targety4)]
    
    print('\n**************************Photometry********************************')
    print('\nPhotometry of Target Star:\n')
    radii = [3., 4., 5.]
    apertures =[CircularAperture(positions, r=r) for r in radii]
    phot_table = aperture_photometry(img_data[0,:,:], apertures)
    target_flux.append(phot_table['aperture_sum_1'][0])
    annulus = [CircularAnnulus(positions, annu_radi_1, annu_radi_2)]
    annu_table = aperture_photometry(img_data[0,:,:], annulus)
    print(annu_table) #we don't need this 
    target_sky.append(annu_table['aperture_sum_0'][0])
    print('Circular Aperture:\n', phot_table)
    print('Circular Annulus:\n', annu_table)
    counts_per_pix = float(annu_table['aperture_sum_0'][0]) / (math.pi * (float(annu_radi_2 ** 2) - float(annu_radi_1 ** 2)))
    print('\nCounts Per Pixel of Taget Star:', counts_per_pix)
    sky_counts = counts_per_pix * (math.pi * (float(radii[0]) ** 2))
    sky_counts_target_lst.append(sky_counts)
    print('\nSky Counts for Target Star:', sky_counts)
    actual_target.append(target_flux[-1] - sky_counts)
    print('\nActual Target:', actual_target[-1])

    print('\nPhotometry of Comparison Star 1:\n')
    apertures_c1 =[CircularAperture(positions_c1, r=r) for r in radii]
    phot_table_c1 = aperture_photometry(img_data[0,:,:], apertures_c1)
    c1_flux.append(phot_table_c1['aperture_sum_1'][0])
    annulus_c1 =[CircularAnnulus(positions_c1, annu_radi_1, annu_radi_2)]
    annu_table_c1 = aperture_photometry(img_data[0,:,:], annulus_c1)
    c1_sky.append(annu_table_c1['aperture_sum_0'][0])
    print('Circular Aperture:\n', phot_table_c1)
    print('Cirular Annulus:\n', annu_table_c1)
    counts_per_pix_1 = float(annu_table_c1['aperture_sum_0'][0]) / (math.pi * (float(annu_radi_2 ** 2) - float(annu_radi_1 ** 2)))
    print('\nCounts Per Pixel of Comparison Star 1:', counts_per_pix_1)
    sky_counts_1 = counts_per_pix_1 * (math.pi * (float(radii[0]) ** 2))
    sky_counts_c1_lst.append(sky_counts_1)
    print('\nSky Counts of Comparison Star 1:', sky_counts_1)
    actual_c1.append(c1_flux[-1] - sky_counts_1)
    print('\nActual Target of Cmparison Star 1:', actual_c1[-1])
    
    print('\nPhotometry of Comparison Star 2:\n')
    apertures_c2 =[CircularAperture(positions_c2, r=r) for r in radii]
    phot_table_c2 = aperture_photometry(img_data[0,:,:], apertures_c2)
    c2_flux.append(phot_table_c2['aperture_sum_1'][0])
    annulus_c2 =[CircularAnnulus(positions_c2, annu_radi_1, annu_radi_2)]
    annu_table_c2 = aperture_photometry(img_data[0,:,:], annulus_c2)
    c2_sky.append(annu_table_c2['aperture_sum_0'][0])
    print('Circular Aperture:\n', phot_table_c2)
    print('Circular Annulus:\n', annu_table_c2)
    counts_per_pix_2 = float(annu_table_c2['aperture_sum_0'][0]) / (math.pi * (float(annu_radi_2 ** 2) - float(annu_radi_1 ** 2)))
    print('\nCounts Per Pixel of Comparison Star 2:', counts_per_pix_2)
    sky_counts_2 = counts_per_pix_2 * (math.pi * (float(radii[0]) ** 2))
    sky_counts_c2_lst.append(sky_counts_2)
    print('\nSky Counts of Comparison Star 2:', sky_counts_2)
    actual_c2.append(c2_flux[-1] - sky_counts_2)
    print('\nActual Target of Cmparison Star 2:', actual_c2[-1])
    
    print('\nPhotometry of Comparison Star 3:\n')
    apertures_c3 =[CircularAperture(positions_c3, r=r) for r in radii]
    phot_table_c3 = aperture_photometry(img_data[0,:,:], apertures_c3)
    c3_flux.append(phot_table_c3['aperture_sum_1'][0])
    annulus_c3 = [CircularAnnulus(positions_c3, annu_radi_1, annu_radi_2)]
    annu_table_c3 = aperture_photometry(img_data[0,:,:], annulus_c3)
    c3_sky.append(annu_table_c3['aperture_sum_0'][0])
    print('Circular Aperture:\n', phot_table_c3)
    print('Circular Annulus:\n', annu_table_c3)
    counts_per_pix_3 = float(annu_table_c3['aperture_sum_0'][0]) / (math.pi * (float(annu_radi_2 ** 2) - float(annu_radi_1 ** 2)))
    print('\nCounts Per Pixel of Comparison Star 3:', counts_per_pix_3)
    sky_counts_3 = counts_per_pix_3 * (math.pi * (float(radii[0]) ** 2))
    sky_counts_c3_lst.append(sky_counts_3)
    print('\nSky Counts of Comparison Star 3:', sky_counts_3)
    actual_c3.append(c3_flux[-1] - sky_counts_3)
    print('\nActual Target of Cmparison Star 3:', actual_c3[-1])
    
    # Now we have to retrieve the reference time in the headers
    print('\n*************************Time Data********************************')
    recorded_time = hdu['OPENTIME']
    print("\nRecorded Time:",recorded_time)
    seperated_times = recorded_time.split(":")
    print('\nTimes sperated:', seperated_times)
    time_in_seconds = float(seperated_times[0]) * 3600 + float(seperated_times[1]) * 60 + float(seperated_times[2])
    print('\nTime in seconds:',time_in_seconds)
    if words == image_list[0]:
        ref_time = time_in_seconds
        time.append(ref_time - ref_time)
        gain = hdu['GAIN']
        rdnoise = hdu['RDNOISE']
    else:
        time.append(time_in_seconds - ref_time)    
    
    print("\nSignal to Noise Ratios:\n")
    print('Target Star:')
    target_v1 = phot_table['aperture_sum_1'][0] * gain
    target_v2 = sky_counts + (rdnoise ** 2) #what is aperture sky counts?
    target_v3 = target_v2 * (math.pi * (radii[0] ** 2))
    target_v4 = target_v1 + target_v3
    target_v5 = target_v4 ** (.5)
    target_sn_ratio = target_v1 / target_v5
    print('Value 1:', target_v1)
    print('Value 2:', target_v2)
    print('value 3:', target_v3)
    print('Value 4:', target_v4)
    print('Value 5:', target_v5)
    print('Signal to Noise Ratio Value of Target Star:', target_sn_ratio)
    
    print('\nComparison Star 1:')
    c1_v1 = phot_table_c1['aperture_sum_1'][0] * gain
    c1_v2 = sky_counts_1 + (rdnoise ** 2)
    c1_v3 = c1_v2 * (math.pi * (radii[0] ** 2))
    c1_v4 = c1_v1 + c1_v3
    c1_v5 = c1_v4 ** (.5)
    c1_sn_ratio = c1_v1 / c1_v5
    print('Value 1:', c1_v1)
    print('Value 2:', c1_v2)
    print('value 3:', c1_v3)
    print('Value 4:', c1_v4)
    print('Value 5:', c1_v5)
    print('Signal to Noise Ratio Value of Comparison Star 1:', c1_sn_ratio)
    
    print('\nComparison Star 2:')
    c2_v1 = phot_table_c2['aperture_sum_1'][0] * gain
    c2_v2 = sky_counts_2 + (rdnoise ** 2)
    c2_v3 = c2_v2 * (math.pi * (radii[0] ** 2))
    c2_v4 = c2_v1 + c2_v3
    c2_v5 = c2_v4 ** (.5)
    c2_sn_ratio = c2_v1 / c2_v5
    print('Value 1:', c2_v1)
    print('Value 2:', c2_v2)
    print('value 3:', c2_v3)
    print('Value 4:', c2_v4)
    print('Value 5:', c2_v5)
    print('Signal to Noise Ratio Value of Comparison Star 1:', c2_sn_ratio)
    
    print('\nComparison Star 3:')
    c3_v1 = phot_table_c1['aperture_sum_1'][0] * gain
    c3_v2 = sky_counts_3 + (rdnoise ** 2)
    c3_v3 = c3_v2 * (math.pi * (radii[0] ** 2))
    c3_v4 = c3_v1 + c3_v3
    c3_v5 = c3_v4 ** (.5)
    c3_sn_ratio = c3_v1 / c3_v5
    print('Value 1:', c3_v1)
    print('Value 2:', c3_v2)
    print('value 3:', c3_v3)
    print('Value 4:', c3_v4)
    print('Value 5:', c3_v5)
    print('Signal to Noise Ratio Value of Comparison Star 1:', c3_sn_ratio)
    
    print('\n*******************************Coordinate Data******************************')
    print('\nX-Axis Coordinates:\n')
    print(time[::-1])

print('\n**********************************Light Curve Draft*******************************')


average_comparison = (np.asarray(c1_flux)+np.asarray(c2_flux)+np.asarray(c3_flux))/3.
#print(target_flux)
#print(average_comparison)
#print(np.median(average_comparison))
#normalized_aperture = (target_flux/average_comparison) / np.median(average_comparison) - 1
normalized_aperture = target_flux/average_comparison



#Remove polynomial
####poly = np.polyfit(time,normalized_aperture,2)
####print(poly)

####print(target_flux)
####print(c1_flux)
####print(normalized_aperture)
####print(time)

a = np.array([time[::-1]])
b = np.array([normalized_aperture])
x = a[0,:]
y = b[0,:]
z = np.polyfit(x, y, 1)
p = np.poly1d(z)
plt.plot(time[::-1], p(time[::-1]), '-')
plt.plot([time[::-1]], [normalized_aperture], 'ro')
plt.title('Light Curve of WD 1150-153',fontsize=16)
plt.ylabel('Fractional Intensity', fontsize=16)
plt.xlabel('Time (Seconds)',fontsize=16)
plt.tick_params(labelsize=16)
#plt.savefig('WD1150-153_lightcurve.png',bbox_inches='tight',dpi=400)

plt.show()

line_division = ([normalized_aperture] / p(time[::-1])) - 1
print('Shape\n', line_division.shape)
print('Shape of P Time', p(time[::-1]).shape)
print('Shpae of Time', time[::-1])
print(line_division)
plt.plot([time[::-1]], line_division, 'bs')
plt.title('Graph of Trendline from Light Curve WD 1150-153\n', fontsize=16)
plt.ylabel("Fractional Instensity of Line Division", fontsize=16)
plt.xlabel("Time (Seconds)", fontsize=16)
plt.show() 

frequency, power = LombScargle(time[::-1], normalized_aperture).autopower()
print(frequency)
print(power)
plt.plot(frequency, power)
plt.ylim(0, 0.05)
plt.xlabel('Frequency', fontsize=16)
plt.ylabel('Power', fontsize=16)
plt.title('Lomb Scargle Chart', fontsize=16)
