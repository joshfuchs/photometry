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

#photutils.test()

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
    hdu = fits.getheader(words) 
    img_data = fits.getdata(words)
    print(img_data.shape)
    targetx, targety = centroid_2dg(img_data[0,348-7:348+7,385-7:385+7])
    print(targetx, targety)
    fig, ax = plt.subplots(1, 1)
    ax.imshow(img_data[0,348-7:348+7,385-7:385+7], origin='lower', interpolation='nearest', cmap='viridis') 
    marker = '+'
    ms, mew = 30, 2.
    plt.plot(targetx, targety, color='#17becf', marker=marker, ms=ms, mew=mew)
    plt.show()
    print(first)
    #componex, componey = centroid_2dg(img_data[56-7:56+7,171-7:171+7])
    #comptwox, comptwoy = centroid_2dg(img_data[129-7:129+7,316-7,316+7])
    #compthreex, compthreey = centroid_2dg(img_data[629-7:629+7,161-7,161+7])









'''
x1, y1 = centroid_2dg(test_image)
x2, y2 = centroid_2dg(test_image)
x3, y3 = centroid_2dg(test_image)
x4, y4 = centroid_2dg(test_image)
fig, ax = plt.subplots('1', '1')
ax.imshow(test_image, origin='lower', interpolation='nearest', cmap='viridis') 
marker = '+'
ms, mew = 30, 2.
plt.plot(x1, y1, color='How?', marker=marker, ms=ms, mew=mew)
plt.plot(x2, y2, color='How?', marker=marker, ms=ms, mew=mew)
plt.plot(x3, y3, color='How?', marker=marker, ms=ms, mew=mew)
plt.plot(x4, y4, color='How?', marker=marker, ms=ms, mew=mew)
# Imports at the beginning


# "Do photometry (photutils)"

print("Step 3:")
   '''








'''
===========================Comments=That=Are=Unneeded==========================
I am proud of myself in this moment
6/22/18
I take it back, I was proud but now I am ashamed 

Dear Monday Morning Roel, you tried your berry best above this line and really
milked whatever is below here. Okay, so you need to incorporate the opening fits
files into the for loop that already prints out the statements of the listphot
file. The output should look like something I don't know. Dear Dr. Fuchs, if you're
reading this, just know that I am trying and if you are reading this, it also
means that pushed this program all by myself :-)
===============================================================================
'''