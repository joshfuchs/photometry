'''
Written by J Meza, updates by J Fuchs. UNC - Chapel Hill

:INPUTS:
       listZero: string, file containing bias images

       listFlat: string, file containing flat images

       listSpec: string, file containing spectral images

       listFe: string, file containing Fe lamp spectra

ReduceSpec.py will automatically sort between Blue and Red ZZ Ceti setups, and will differentiate between different targets. But each observation of a target must be grouped in the file.

:OUTPUTS:
       File structure:
         - b*fits: bias-subtracted images
         - fb*fits: images flat-fielded after bias subtraction
         - tfb*fits: images trimmed after flat-fielding and bias subtraction

       reduction_DATE.txt: file containing info for diagnostic purposes. DATE is date and time of reduction. Zeros in a whole column typically mean blue or red setup not included. Columns are: 0) average from bias, 1) average from scaled bias, 2) standard deviation of bias, 3) Blue flat field average, 4) Blue flat field standard deviation, 5) Blue flat field scaled average, 6) Blue flat field scaled standard deviation, 7) Red flat field average, 8) Red flat field standard deviation, 9) Red flat field scaled average, 10) Red flat field scaled standard deviation, 11) Combined blue flat pixel values, 12) Combined blue flat values, 13) Polynomial fit to combined blue flat \n 14) Combined red flat pixel values, 15) Combined red flat values, 16) Polynomial fit to combined red flat

       FWHM_records_DATE.txt: file containing information about the FWHM and profile locations for 2D images in listSpec. Columns are: filename, Column of 2D image checked, FWHM of Gaussian fit to that column, Center position of Gaussian fit to that column, Second column checked, FWHM of second column, Center position of second column.


'''

# ===========================================================================
# Packages ==================================================================
# ===========================================================================

import numpy as np
import os
import ReduceSpec_tools as rt
import warnings

# ===========================================================================
# Code to Reduce Spectrum ===================================================
# ===========================================================================
# This peice takes in arguments from the command line to 
# All the functions required are called from ReduceSpec_tools.py



def reduce_now(args):
    nargs = len(args)
    if (nargs < 4):
        print("\n====================\n")
        print("\nNot Enough Inputs." )
        print("Need at least 4 inputs: listZero, listFlat, listSpec, listFe")
        print("Optional inputs: overwrite= , low_sig= , high_sig=  ")
        print("Example:")
        print("\n>>> python imcombine.py listZero listFlat listSpec listFe \n")
        print("\n====================\n")
    
    # Unpack list from command line and combe trough them for diffrent observations # 
    scriptname = args[0]
    zero_lists = rt.List_Combe( rt.Read_List( args[1] ) )
    flat_lists = rt.List_Combe( rt.Read_List( args[2] ) )
    spec_lists = rt.List_Combe( rt.Read_List( args[3] ) )
    
    # Select names from the first image of each observation # 
    zero_names= []
    for zero in zero_lists:
        zero_names.append(zero[0][5:])
    flat_names= []
    for flat in flat_lists:
        flat_names.append(flat[0][5:])
    spec_names= []
    for spec in spec_lists:
        spec_names.append(spec[0][5:])
    
    # Default values for special commands if none are given these dont change #   
    overwrite = True # dont give imcombine permision to overwrite files # 
    lo_sig = 10
    hi_sig = 3
    method = 'median' # method used to combine images 
    
    # If overwrite special comand is given # 
    if nargs >= 6:
        overwrite = args[5]
        warnings.filterwarnings('ignore', category=UserWarning, append=True)
    # If low_sigma and high_sigma values are given # 
    if nargs >= 8: 
        lo_sig = float(args[6])
        hi_sig = float(args[7]) 
    # If method is given #  
    if nargs >= 9:
        method = args[8]
        
    #Set up array to save for diagnostics. This is defined in rt.init()
    rt.init()


    
    # The rest of the code runs the reduction procces up to apall #  =========
    # Combine Zeros # 
    comb_zero = rt.imcombine(zero_lists[0], zero_names[0], 'average', lo_sig= 10, 
                        hi_sig= 3, overwrite= overwrite)
    
    # Bias Subtract Flats # 
    nf= len(flat_lists) # number of flats
    b_flat_lists= []
    i= 0
    while i < nf:
        b_flat_lists.append( rt.Bias_Subtract(flat_lists[i], comb_zero ) )
        i= i+1
    
    # Combine Bias Subtracted Flats # 
    i= 0
    comb_flat= []
    while i < nf:
        comb_flat.append( rt.imcombine(b_flat_lists[i], 'b.'+flat_names[i], 'median', 
                    lo_sig= 10, hi_sig= 3, overwrite= overwrite) )
        i= i+1
    
    
    # Normalize Flat # 
    i= 0
    nb_flat= []
    while i < nf:
        nb_flat.append( rt.Norm_Flat_Avg(comb_flat[i])) # (divide by average of counts)
        i= i+1
  
    
    print('tennisten')
    # Bias Subtract Spec # 
    i= 0
    b_spec_list= []
    nsp= len(spec_lists); # number of spectra
    while i < nsp:
        b_spec_list.append( rt.Bias_Subtract(spec_lists[i], comb_zero) )
        i= i+1
        
        
    # Flat Field Individual Spectra #
    i= 0
    ftb_spec_list = []

    tb_spec_list = rt.List_Combe(b_spec_list)
  
    print(tb_spec_list[i])
    print(type(tb_spec_list[i]))

    #tb_spec_list = rt.List_Combe(b_spec_list)
    tb_spec_list = b_spec_list

    
    while i < nsp:
        ftb_spec_list.append( rt.Flat_Field(tb_spec_list[i], nb_flat[0]) )
        i = i +1

    '''
    blueindex = [i for i, s in enumerate(nb_flat) if 'blue' in s.lower()]
    nbflatblue = nb_flat[blueindex[0]]
    redindex = [i for i, s in enumerate(nb_flat) if 'red' in s.lower()]
    if len(redindex) > 0:
        nbflatred = nb_flat[redindex[0]]
    i= 0
    ftb_spec_list = []
    tb_spec_list = rt.List_Combe(b_spec_list)
    while i < nsp:
        if tb_spec_list[i][0].lower().__contains__('blue') == True:
            ftb_spec_list.append( rt.Flat_Field(tb_spec_list[i], nbflatblue) )
        elif tb_spec_list[i][0].lower().__contains__('red') == True:
            ftb_spec_list.append( rt.Flat_Field(tb_spec_list[i], nbflatred) )
        else: 
            print("Problem applying the Flats.")
            print("Could not identify blue or red setup.")
        i= i+1
         tb_spec_list = rt.List_Combe(b_spec_list)
  
    print(tb_spec_list[i])
    print(type(tb_spec_list[i]))
    '''
    
    
    # Save all diagnostic info
    rt.save_diagnostic()
    
    #LA Cosmic
    i = 0
    cftb_spec = []
    cftb_mask = []
    while i < nsp:
        m = 0
        while m < len(ftb_spec_list[i]):
            lacos_spec, lacos_mask = rt.lacosmic(ftb_spec_list[i][m])
            cftb_spec.append(lacos_spec)
            cftb_mask.append(lacos_mask)
            m += 1
        i += 1
    
    cftb_spec_list = rt.List_Combe(cftb_spec)
    cftb_mask_list = rt.List_Combe(cftb_mask)
    
   
   

    print("Done. Ready for Apeture Extraction.\n")
    
# ===========================================================================
# ===========================================================================
# ===========================================================================

#To run from command line
if __name__ == "__main__":
    from sys import argv
    #args = argv # arguments from comand line #
    args = ['script','listzero','listflat','listphot']
    #os.chdir('C:/Users/rrolvera/Astro Research/data/WD1150-153/photometry')
    reduce_now(args)
