# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 11:18:53 2018

@author: rrolvera
"""

import numpy as np # You always have to do this, the np is to avoid confusion
import os # You also have to do this, it stand for operating system
from astropy.io import fits # Imported so we can open the fits files
from photutils import centroid_2dg, CircularAperture, SkyCircularAperture # for the centroid
import matplotlib.pyplot as plt # Imports for the plots
from astropy import units as u 
from astropy.coordinates import SkyCoord
from photutils import aperture_photometry, make_source_mask
from astropy.stats import sigma_clipped_stats
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize

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

# Now we can write a "for" loop 

print("Part 2:")

# First, we want to write a "for" loop that reads out the list as statements
# Read in Fits File (astropy.io.fits)
# How in the universe to do insert open.fits into this loop?

print("Step 1:")

images = image_list # This defines that the word dummy is the dummy_lists
for words in images: # This "for" loop is for the statements in dummy_lists
    hdu = fits.getheader(words) # This gets the words from the header
    img_data = fits.getdata(words) # This reads it as data
    print('Image Dimensions:')
    print(img_data.shape) # This shows us the the list
    print('Diagrams:')
    marker = '+'
    ms, mew = 30, 2.
    targetx1, targety1 = centroid_2dg(img_data[0,y1-7:y1+7,x1-7:x1+7])
    print('Centroid Target Star:')
    fig, ax = plt.subplots(1, 1)
    ax.imshow(img_data[0,y1-7:y1+7,x1-7:x1+7], origin='lower', interpolation='nearest', cmap='viridis')
    plt.plot(targetx1, targety1, color='#17becf', marker=marker, ms=ms, mew=mew)
    plt.show()
    print(targetx1, targety1)
    targetx2, targety2 = centroid_2dg(img_data[0,y2-7:y2+7,x2-7:x2+7])
    print('Centroid Comparison Star 1:')
    fig, ax = plt.subplots(1, 1)
    ax.imshow(img_data[0,y2-7:y2+7,x2-7:x2+7], origin='lower', interpolation='nearest', cmap='viridis')
    plt.plot(targetx2, targety2, color='r', marker=marker, ms=ms, mew=mew)
    plt.show()
    print(targetx2, targety2)
    targetx3, targety3 = centroid_2dg(img_data[0,y3-7:y3+7,x3-7:x3+7])
    print('Centroid Comparison Star 2:')
    fig, ax = plt.subplots(1, 1)
    ax.imshow(img_data[0,y3-7:y3+7,x3-7:x3+7], origin='lower', interpolation='nearest', cmap='viridis')
    plt.plot(targetx3, targety3, color='y', marker=marker, ms=ms, mew=mew)
    plt.show()
    print(targetx3, targety3)
    targetx4, targety4 = centroid_2dg(img_data[0,y4-7:y4+7,x4-7:x4+7])
    print('Centroid Comparison Star 3:')
    fig, ax = plt.subplots(1, 1)
    ax.imshow(img_data[0,y4-7:y4+7,x4-7:x4+7], origin='lower', interpolation='nearest', cmap='viridis')
    plt.plot(targetx4, targety4, color='k', marker=marker, ms=ms, mew=mew)
    plt.show()
    print(targetx4, targety4)
    positions = [(30., 30.), (40., 40.)]
    apertures = CircularAperture(positions, r=3)
    positions = SkyCoord(l=[1.2, 2.3] * u.deg, b=[0.1, 0.2] * u.deg, frame='galactic')
    apertures = SkyCircularAperture(positions, r=4. * u.arcsec)
    positions = [(x1 + targetx1, y1 + targety1), (x1 + targetx1, y1 + targety1)]
    apertures = CircularAperture(positions, r=3)
    data = np.ones((100, 100))
    mean, median, std = sigma_clipped_stats(data, sigma=3.0, iters=5)
    norm = ImageNormalize(stretch=SqrtStretch())
    mask = make_source_mask(data, snr=2, npixels=5, dilate_size=11)
    mean, median, std = sigma_clipped_stats(data, sigma=3.0, mask=mask)
    print(mean, median, std)
    ny, nx = data.shape
    y, x = np.mgrid[:ny, :nx]
    gradient = x * y / 5000.
    data2 = data + gradient
    phot_table = aperture_photometry(data, apertures, method='subpixel', subpixels=5)
    # Multiple Apertures at Each Position
    print('Multiple Apertures at Each Position:')
    radii = [3., 4., 5.]
    apertures =[CircularAperture(positions, r=r) for r in radii]
    phot_table = aperture_photometry(img_data[0,:,:], apertures)
    print(phot_table)
    print('Isolated Value:')
    print(phot_table['aperture_sum_1'][0])
    # Now we can add this value to the "target_flux" we initialized earlier
    print('Y-Axis Coordinates:')
    target_flux.append(phot_table['aperture_sum_1'][1])
    print(target_flux)
    # Now we have to retrieve the reference time in the headers
    recorded_time = hdu['OPENTIME']
    print("Recorded Time")
    print(recorded_time)
    seperated_times = recorded_time.split(":")
    print(seperated_times)
    time_in_seconds = float(seperated_times[0]) * 3600 + float(seperated_times[1]) * 60 + float(seperated_times[2])
    print(time_in_seconds)
    if words == image_list:
        ref_time = time_in_seconds
        time.append(ref_time - ref_time)
    else:
        time.append(ref_time - time_in_seconds)
    gain = hdu['GAIN']
    rdnoise = hdu['RDNOISE']
    print("Gain:")
    print(gain)
    print("RD Noise")
    print(rdnoise)
    print("Signal to Noise Ratio:")
    value_1 = phot_table['aperture_sum_1'][0] * gain
    print("Value 1:")
    print(value_1)
    print(mean, median, std)
print('Y-Axis Coordinates:')
target_flux.append(phot_table['aperture_sum_1'][1])
print(target_flux)
print('X-Axis Coordinates:')
print(time[::-1])    

print('Part Three:')

print('Light Curve Draft:')

normalized_aperture = target_flux / np.median(target_flux) - 1

plt.plot([time[::-1]], [normalized_aperture], 'ro')
plt.title('Light Curve of Target Star')
plt.ylabel('Fractional Intensity')
plt.xlabel('Time (Seconds)')
plt.show()