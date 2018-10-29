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

os.chdir('C:/Users/rrolvera/Astro Research/data/WD1150-153/photometry') # Changes our directory to where our "listphot" is located
images = 'listphot' # A name is given to the file we are using/want 
image_list = Read_List(images) # Read_List command reads dummy and is called dummy_lists as a whole
print(image_list) # This reads out our list

# Second, we want to define the four stars we will be using
# "Set coordinates of initial target + 3 comparisons"

print("Step 2:")

target_star = [] # This lets us set a name to a set of coordinates
x1 = 385 # This is the x coordinate of the star 
y1 = 348 # This is the y coordinate of the star
target_star.append((x1, y1)) # This adds our coordinate/star to our program
print("Target Star Coordinates:") # This is the title of information
print(target_star) # This prints our star coordinates

comparison_star_one = [] # This lets us set a name to a set of coordinates
x2 = 56 # This is the x coordinate of the star 
y2 = 171 # This is the y coordinate of the star
comparison_star_one.append((x2, y2)) # This adds our coordinate/star to our program
print("Comparison Star One Coordinates:") # This is the title of information
print(comparison_star_one) # This prints our star coordinates
        
comparison_star_two = [] # This lets us set a name to a set of coordinates
x3 = 129 # This is the x coordinate of the star 
y3 = 316 # This is the y coordinate of the star
comparison_star_two.append((x3, y3)) # This adds our coordinate/star to our program
print("Comparison Star Two Coordinates:") # This is the title of information
print(comparison_star_two) # This prints our star coordinates

comparison_star_three = [] # This lets us set a name to a set of coordinates
x4 = 629 # This is the x coordinate of the star 
y4 = 161 # This is the y coordinate of the star
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
# Now we can write a "for" loop 

print("Part 2:")

# First, we want to write a "for" loop that reads out the list as statements
# Read in Fits File (astropy.io.fits)
# How in the universe to do insert open.fits into this loop?

print("Step 1:")
pix_dis = 7
images = image_list # This defines that the word dummy is the dummy_lists
for words in images: # This "for" loop is for the statements in dummy_lists
    hdu = fits.getheader(words) # This gets the words from the header
    img_data = fits.getdata(words) # This reads it as data
    print('Image Dimensions:')
    print(img_data.shape) # This shows us the the list
    print('Diagrams:')
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
    print('Centroid Comparison Star 1:')
    fig, ax = plt.subplots(1, 1)
    ax.imshow(img_data[0,y2-pix_dis:y2+pix_dis,x2-pix_dis:x2+pix_dis], origin='lower', interpolation='nearest', cmap='viridis')
    plt.plot(targetx2, targety2, color='r', marker=marker, ms=ms, mew=mew)
    plt.show()
    print(targetx2, targety2)
    targetx3, targety3 = centroid_2dg(img_data[0,y3-pix_dis:y3+pix_dis,x3-pix_dis:x3+pix_dis])
    print('Centroid Comparison Star 2:')
    fig, ax = plt.subplots(1, 1)
    ax.imshow(img_data[0,y3-pix_dis:y3+pix_dis,x3-pix_dis:x3+pix_dis], origin='lower', interpolation='nearest', cmap='viridis')
    plt.plot(targetx3, targety3, color='y', marker=marker, ms=ms, mew=mew)
    plt.show()
    print(targetx3, targety3)
    targetx4, targety4 = centroid_2dg(img_data[0,y4-pix_dis:y4+pix_dis,x4-pix_dis:x4+pix_dis])
    print('Centroid Comparison Star 3:')
    fig, ax = plt.subplots(1, 1)
    ax.imshow(img_data[0,y4-pix_dis:y4+pix_dis,x4-pix_dis:x4+pix_dis], origin='lower', interpolation='nearest', cmap='viridis')
    plt.plot(targetx4, targety4, color='k', marker=marker, ms=ms, mew=mew)
    plt.show()
    print(targetx4, targety4)
    positions = [(x1 + targetx1, y1 + targety1), (x1 + targetx1, y1 + targety1)]
    positions_c1 =  [(x2 + targetx2, y2 + targety2)]
    positions_c2 = [(x3 + targetx3, y3 + targety3)]
    positions_c3 = [(x4 + targetx4, y4 + targety4)]
    # Multiple Apertures at Each Position
    print('Multiple Apertures at Each Position:')
    radii = [3., 4., 5.]
    apertures =[CircularAperture(positions, r=r) for r in radii]
    phot_table = aperture_photometry(img_data[0,:,:], apertures)
    print(phot_table)
    print('Circular Annulus:')
    annulus = [CircularAnnulus(positions, 8, 11)]
    annu_table = aperture_photometry(img_data[0,:,:], annulus)
    print(annu_table)
    print('Isolated Value:')
    print(phot_table['aperture_sum_1'][0])
    # Now we can add this value to the "target_flux" we initialized earlier
    print('Y-Axis Coordinates:')
    target_flux.append(phot_table['aperture_sum_1'][0])
    print(target_flux)
    #Do photometry of comparison stars
    print('Photometry of Comparison Star 1: ')
    apertures_c1 =[CircularAperture(positions_c1, r=r) for r in radii]
    phot_table_c1 = aperture_photometry(img_data[0,:,:], apertures_c1)
    c1_flux.append(phot_table_c1['aperture_sum_1'][0])
    annulus_c1 =[CircularAnnulus(positions_c1, 8, 11)]
    annu_table_c1 = aperture_photometry(img_data[0,:,:], annulus)
    print('Circular Aperture: ', phot_table_c1)
    print('Cirular Annulus:\n', annu_table_c1)
    
    print('Photometry of Comparison Star 2: ')
    apertures_c2 =[CircularAperture(positions_c2, r=r) for r in radii]
    phot_table_c2 = aperture_photometry(img_data[0,:,:], apertures_c2)
    c2_flux.append(phot_table_c2['aperture_sum_1'][0])
    annulus_c2 =[CircularAnnulus(positions_c2, 8, 11)]
    annu_table_c2 = aperture_photometry(img_data[0,:,:], annulus)
    print('Circular Aperture: ', phot_table_c2)
    print('Circular Annulus:\n', annu_table_c2)
    
    print('Photometry of Comparison Star 3 :')
    apertures_c3 =[CircularAperture(positions_c3, r=r) for r in radii]
    phot_table_c3 = aperture_photometry(img_data[0,:,:], apertures_c3)
    c3_flux.append(phot_table_c3['aperture_sum_1'][0])
    annulus_c3 = [CircularAnnulus(positions_c3, 8, 11)]
    annu_table_c3 = aperture_photometry(img_data[0,:,:], annulus)
    print('Circular Aperture: ', phot_table_c3)
    print('Circular Annulus: ', annu_table_c3)
    
    # Now we have to retrieve the reference time in the headers
    recorded_time = hdu['OPENTIME']
    print("Recorded Time",recorded_time)
    seperated_times = recorded_time.split(":")
    print(seperated_times)
    time_in_seconds = float(seperated_times[0]) * 3600 + float(seperated_times[1]) * 60 + float(seperated_times[2])
    print('Time in seconds:',time_in_seconds)
    if words == image_list[0]:
        ref_time = time_in_seconds
        time.append(ref_time - ref_time)
        gain = hdu['GAIN']
        rdnoise = hdu['RDNOISE']
    else:
        time.append(time_in_seconds - ref_time)    
    
    print("Signal to Noise Ratio:")
    value1 = phot_table['aperture_sum_1'][0] * gain
    value2 = phot_table['aperture_sum_1'][0] + (rdnoise ** 2) #what is aperture sky counts?
    value3 = value2 * (math.pi * (radii[0] ** 2))
    value4 = value1 + value3
    value5 = value4 ** (.5)
    sn_ratio = value1 / value5
    print('Value 1:',value1)
    print('Value 2:',value2)
    print('value 3:',value3)
    print('Value 4:',value4)
    print('Value 5:',value5)
    print('Ratio:',sn_ratio)
    print(target_flux)
    print('X-Axis Coordinates:')
    print(time[::-1])    

print('Part Three:')

print('Light Curve Draft:')


average_comparison = (np.asarray(c1_flux)+np.asarray(c2_flux)+np.asarray(c3_flux))/3.
#print(target_flux)
#print(average_comparison)
#print(np.median(average_comparison))
#normalized_aperture = (target_flux/average_comparison) / np.median(average_comparison) - 1
normalized_aperture = (target_flux/average_comparison) / np.median((target_flux/average_comparison)) - 1

#Remove polynomial
####poly = np.polyfit(time,normalized_aperture,2)
####print(poly)

####print(target_flux)
####print(c1_flux)
####print(normalized_aperture)
####print(time)
plt.plot([time[::-1]], [normalized_aperture], 'ro')
plt.title('Light Curve of WD 1150-153',fontsize=16)
plt.ylabel('Fractional Intensity', fontsize=16)
plt.xlabel('Time (Seconds)',fontsize=16)
plt.tick_params(labelsize=16)
#plt.savefig('WD1150-153_lightcurve.png',bbox_inches='tight',dpi=400)

plt.show()