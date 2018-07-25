# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 11:18:53 2018

@author: rrolvera
"""

import numpy as np # You always have to do this, the np is to avoid confusion
import os # You also have to do this, it stand for operating system
from astropy.io import fits # Imported so we can open the fits files
from photutils import centroid_com, centroid_1dg, centroid_2dg # for the centroid
import matplotlib.pyplot as plt # Imports for the plots
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes # Tool
from mpl_toolkits.axes_grid1.inset_locator import mark_inset # Tool
import matplotlib.pyplot 
from photutils import CircularAperture
from astropy import units as u 
from astropy.coordinates import SkyCoord
from photutils import SkyCircularAperture 
from photutils import aperture_photometry 
from astropy.stats import SigmaClip
from photutils import Background2D, MedianBackground
from photutils import CircularAnnulus
from astropy.stats import mad_std
from astropy.stats import sigma_clipped_stats
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.stats import biweight_location
from photutils import make_source_mask
from photutils import MeanBackground


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
dummy = 'listphot' # A name is given to the file we are using/want 
dummy_lists = Read_List( dummy ) # Read_List command reads dummy and is called dummy_lists as a whole
print(dummy_lists) # This reads out our list

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
# "Now we can initiliaze arrays: time, target flux, comp flux, errors: length of number of files"
# Did I even do this the right way? Dallas and Andrew said this is how you initialize
print("Step 3:")

time = [] # This initilizes the array, in the future we can input information
target_flux = [] # This initilizes the array, in the future we can input information
comp_flux = [] # This initilizes the array, in the future we can input information
errors = [] # This initilizes the array, in the future we can input information
file_amount = [] # This initilizes the array, in the future we can input information

# Now we can write a "for" loop 

print("Part 2:")

# First, we want to write a "for" loop that reads out the list as statements
# Read in Fits File (astropy.io.fits)
# How in the universe to do insert open.fits into this loop?

print("Step 1:")

dummy = dummy_lists # This defines that the word dummy is the dummy_lists
for words in dummy: # This "for" loop is for the statements in dummy_lists
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
    # Creating Apeture Objects
    #print('Creating Aperture Objects:')
    positions = [(30., 30.), (40., 40.)]
    apertures = CircularAperture(positions, r=3)
    positions = SkyCoord(l=[1.2, 2.3] * u.deg, b=[0.1, 0.2] * u.deg, frame='galactic')
    apertures = SkyCircularAperture(positions, r=4. * u.arcsec)
    # Performing Aperture Photometry
    #print('Performing Aperture Photometry:')
    positions = [(30., 30.), (40., 40.)]
    apertures = CircularAperture(positions, r=3)
    data = np.ones((100, 100))
    mean, median, std = sigma_clipped_stats(data, sigma=3.0, iters=5)
    norm = ImageNormalize(stretch=SqrtStretch())
    #plt.imshow(data, norm=norm, origin='lower', cmap='Greys_r')
    #print(np.median(data))
    #print(biweight_location((data)))
    #print(mad_std)
    #print(mean, median, std)
    mask = make_source_mask(data, snr=2, npixels=5, dilate_size=11)
    mean, median, std = sigma_clipped_stats(data, sigma=3.0, mask=mask)
    print(mean, median, std)
    ny, nx = data.shape
    y, x = np.mgrid[:ny, :nx]
    gradient = x * y / 5000.
    data2 = data + gradient
    #plt.imshow(data2, norm=norm, origin='lower', cmap='Greys_r')
    #sigma_clip = 
    #phot_table = aperture_photometry(data, apertures)
    #print(phot_table)
    # Apeture and Pixel Overlap
    #print('Aperture and Pixel Overlap:')
    phot_table = aperture_photometry(data, apertures, method='subpixel', subpixels=5)
    #print(phot_table)
    # Multiple Apertures at Each Position
    print('Multiple Apertures at Each Position:')
    radii = [3., 4., 5.]
    apertures =[CircularAperture(positions, r=r) for r in radii]
    phot_table = aperture_photometry(data, apertures)
    print(phot_table)
    # Before we go on, we have to establish the backgroun
    sigma_clip = SigmaClip(sigma=3.0, iters=10)
    bkg_estimator = MedianBackground(sigma_clip)
    #bkg = Background2D(data, (50, 50), filter_size=(3, 3), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
    bkg = MeanBackground(sigma_clip)
    bkg_value = bkg(data)
    #print(bkg_value)
    #print(bkg)
    #print(bkg_estimator)
    #print(data)
    #print(apertures)
    # Global Background Subtractio
    #print('Global Background Subtraction:')
    #phot_table = aperture_photometry(data - bkg, apertures) # Can be ran without "data - bkg, apertures"
    #print(phot_table)
    # Local Background Subtraction
    #print('Local Background Subtraction:')
    apertures = CircularAperture(positions, r=3)
    annulus_apertures = CircularAnnulus(positions, r_in=6., r_out=8)
    # Perform photometry in both apertures
    #print('Perfrom Photometry in Both Apertures:')
    apers = [apertures, annulus_apertures]
    phot_table = aperture_photometry(data, apers)
   # print(phot_table)
    # Calculate the mean local background within the circular annulus aperture
    #print('Calculate the Mean Local Background Within the Circular Annulus Aperture:')
    bkg_mean = phot_table['aperture_sum_1'] / annulus_apertures.area()
    bkg_sum = bkg_mean * apertures.area()
    final_sum = phot_table['aperture_sum_0'] - bkg_sum
    phot_table['residual_aperture_sum'] = final_sum
    #print(phot_table['residual_aperture_sum'])
    # Error Estimation
   # print('Error Estimation:')
    error = 0.1 * data
    phot_table = aperture_photometry(data, apertures, error=error)
   # print(phot_table)
    # Pixel Masking
    #print('Pixel Masking:')
    data = np.ones((5, 5))
    aperture = CircularAperture((2, 2), 2.)
    mask = np.zeros_like(data, dtype=bool)
    data[2, 2] = 100. # Bad pixel
    mask[2, 2] = True
    t1 = aperture_photometry(data, aperture, mask=mask)
    print(t1['aperture_sum'])
    # Isolate the value 
    print('Isolated Value:')
    print(phot_table['aperture_sum'][0])
    # Now we can add this value to the "target_flux" we initialized earlier
    print('Y-Axis Coordinates:')
    target_flux.append(phot_table['aperture_sum'][1])
    print(target_flux)
    # Now we have to retrieve the reference time in the headers
    spec_data = fits.getdata(words)
    spec_header = fits.getheader(words)
    raw_time = spec_header['OPENTIME']
    print("Raw Time")
    print(raw_time)
    diff_times = raw_time.split(":")
    print(diff_times)
    time1 = float(diff_times[0]) * 3600 + float(diff_times[1]) * 60 + float(diff_times[2])
    print(time1)
    if words == dummy_lists:
        ref_time = time1
        time.append(ref_time - ref_time)
    else:
        time.append(ref_time - time1)
    gain = spec_header['GAIN']
    rdnoise = spec_header['RDNOISE']
    print("Gain:")
    print(gain)
    print("RD Noise")
    print(rdnoise)
    print("Signal to Noise Ratio:")
    value_1 = phot_table['aperture_sum'][0] * gain
    print("Value 1:")
    print(value_1)
    print(mean, median, std)
    #value_2 = 
print('Y-Axis Coordinates:')
target_flux.append(phot_table['aperture_sum'][1])
print(target_flux)
print('X-Axis Coordinates:')
print(time[::-1])    

print('Part Three:')

print('Light Curve Draft:')

plt.plot([time[::-1]], [target_flux], 'ro')
plt.title('Light Curve of Target Star')
plt.ylabel('Luminosity')
plt.xlabel('Time')
plt.show()

end = ("Progress")
for words in end:
    print(words)
    