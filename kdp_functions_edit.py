#!/usr/bin/env python2.7

## @file
#  @brief   kdp_functions.py contains functions required for single scan kdp generation
#  @program kdp_functions.py contains functions required for single scan kdp generation

# - S. Best, September '15

# -----------------------------------------------------------------------------
# Imports
# -----------------------------------------------------------------------------

import os
import sys
import numpy as np
import scipy.signal as sg

# libradarnet4 directory (this can be removed):
#sys.path.append(os.getenv('RADARNET4_BASE') + '/lib')   
#import libradarnet4 as rn
    
# ----------------------------------------------------------------------------
# extract quantities:
# ----------------------------------------------------------------------------

def extract_data(inFile):
    """Extract dual pol quantities + ray / bin info."""
    
    # create polar volume from input file:
    pv_in = rn.PolarVolume(inFile)
    
    # standard fields:
    zv = pv_in.getData('ZV')                          # reflectivity (linear units)
    rays = pv_in.getScanNRays(0)                      # number of rays
    bins = pv_in.getScanNBins(0)                      # number of bins 
    binLength = pv_in.getScanBinLengthKM(0) * 1000.0  # bin length (metres)
    flags = pv_in.getData('FLAGS')                    # QC flags
    
    # get dual pol parameters (exit if not found):
    try:
        phidp = pv_in.getData('PHIDP')
    except rn.RadarnetError as detail:
        print '\n*** INFO: exception occured -', detail
        print ' - **EXITING**'
        sys.exit(0) # this is not an error condition because we may have singlepol input        
    try:
        rhohv = pv_in.getData('RHOHV')
    except rn.RadarnetError as detail:
        print '\n*** INFO: exception occured -', detail
        print ' - **EXITING**'
        sys.exit(0) # this is not an error condition because we may have singlepol input       
    try:    
        zdr = pv_in.getData('ZDR')
    except rn.RadarnetError as detail:
        print '\n*** INFO: exception occured -', detail
        print ' - **EXITING**'
        sys.exit(0) # this is not an error condition because we may have singlepol input
                 
    return pv_in, zv, phidp, rhohv, zdr, rays, bins, binLength, flags
           
# -----------------------------------------------------------------------------
# generate meteo mask:
# -----------------------------------------------------------------------------

def generate_meteo_mask(rays, bins, flags, rhohv, meteoThresh):
    """Meteo mask is used for intial clutter removal from phidp."""
    
    # initially set all to one:
    meteoMask = np.ones((rays, bins))
    
    # set to zero where below rho_hv meteo threshold:
    meteoMask[np.where(rhohv < meteoThresh)] = 0
    
    # set mask to zero where flags are non-zero (i.e. not valid weather):
    meteoMask[(flags != 0)] = 0
    
    # if no valid pixels in entire scan, exit:
    validPix = meteoMask[(meteoMask==1)]
    if np.size(validPix) == 0:
        print '\n*** INFO: No meteo pixels detected - skipping KDP calculation for this scan. ***\n'
        sys.exit(0) # this is not an error condition
        
    return meteoMask
    
# -----------------------------------------------------------------------------
# generate rain mask:
# -----------------------------------------------------------------------------
    
def generate_rain_mask(rays, bins, rhohv, rainThresh):
    """Rain mask is to remove non-rain regions before calculating KDP."""
    
    # initially set all to one:
    rainMask = np.ones((rays, bins))
    
    # set mask to zero where rhohv is below rain threshold:
    rainMask[(rhohv < rainThresh)] = 0
    
    # if no valid pixels in entire scan, exit:
    validPix = rainMask[(rainMask==1)]
    if np.size(validPix) == 0:
        print '\n*** INFO: No rain pixels detected - skipping KDP calculation for this scan. ***\n'
        sys.exit(0) # this is not an error condition
    
    return rainMask
    
# -----------------------------------------------------------------------------
# unwrap phidp:
# ----------------------------------------------------------------------------- 

def unwrap_phidp(rays, bins, mask, phidp):
    """Carries out phase-unwrapping of phidp
    
    (Phase-wrapping may affect phidp under certain circumstances / 
    with certain radars).
    """

    # initially set to copy of existing array:
    phidp_unwrap = np.copy(phidp)

    # loop over rays:
    for i in range(0, rays):
    
        # only continue for rays with clean bins:
        cleanBins = mask[i][(mask[i]==1)]
        if len(cleanBins) > 0:
        
            # set last (most recent) clean pixel along ray:
            j_lastClean = 0
        
            # loop over bins:
            for j in range(1, bins):
            
                # only look for clean pixels:
                if (mask[i,j] == 1):
                
                    # if apparent change is big enough:
                    if abs(phidp_unwrap[i,j] - phidp_unwrap[i,j_lastClean]) \
                        > 340.0:
                        
                        # calculate 'actual' change:
                        delta = (180.0 - abs(phidp_unwrap[i,j_lastClean])) + \
                                (180.0 - abs(phidp_unwrap[i,j]))
                    
                        # if first value is negative, make second one negative:
                        if np.sign(phidp_unwrap[i,j_lastClean]) == -1:
                            phidp_unwrap[i,j] = phidp_unwrap[i,j_lastClean] \
                            - delta
                        # if first value is positive, make second one positive:
                        else:
                            phidp_unwrap[i,j] = phidp_unwrap[i,j_lastClean] \
                            + delta
                            
                    # reset latest clean pixel:        
                    j_lastClean = j
                        
    return phidp_unwrap
                   
# -----------------------------------------------------------------------------
# clean phidp:
# -----------------------------------------------------------------------------    
    
def clean_phidp(rays, bins, phidp, mask, filterBins):
    """Removes flagged values from phidp, rebuilds profile using interpolation.
    
    (Continuous phid_dp is required for gradient calculation).
    """
    
    # bin range:
    binRange = np.arange(0, bins)

    # pre-define cleaned phidp: 
    phidp_clean = np.zeros((rays, bins))

    # loop over rays:
    for i in range(0, rays):  
        
        # only continue for rays with clean bins:
        cleanBins = mask[i][(mask[i]==1)]
        if len(cleanBins) > 0:
        
            # exctract phidp values where unmasked:        
            phidp_unmasked = phidp[i][(mask[i] == 1)]
            
            # determine filter length 
            # (use pre-defined value where possible):
            if len(cleanBins) >= filterBins:
                fb = filterBins
            else: fb = ((len(cleanBins)) // 2) + 1
            # amend if even:
            if fb % 2 == 0:
                fb += 1
                
            # do filtering:
            phidp_unmasked_filt = sg.medfilt(phidp_unmasked, fb)
            # compensate for edge effects:
            phidp_unmasked_filt[0:fb//2] = phidp_unmasked_filt[fb//2]
            phidp_unmasked_filt[-fb//2: ] = phidp_unmasked_filt[-(fb//2)-1]
            # assign to final array:
            phidp_clean[i][(mask[i] == 1)] = phidp_unmasked_filt
            # set up interpolation vectors for current ray:
            x_clean = binRange[(mask[i] == 1)]
            y_clean = phidp_clean[i][(mask[i] == 1)]
            # do interpolation:
            phidp_clean[i][(mask[i] == 0)] = \
            np.interp(binRange[(mask[i] == 0)], x_clean, y_clean)
                                                                                                      
    return phidp_clean      
    
# -----------------------------------------------------------------------------
# calculate kdp (method 2 - faster):
# -----------------------------------------------------------------------------

def calc_kdp_v2(rays, bins, binLength, phidp, mask):
    """Calculates KDP as the gradient of reconstucted phidp along each ray
    
    Uses numpy.gradient function - quicker than linear regression

    >>> import numpy as np; calc_kdp_v2(2,2,.6,np.array([[0,6],[0,6]]), np.ones((2,2)))
    info: ignoring rays with no clean pixels ...
    array([[ 5000.,  5000.],
           [ 5000.,  5000.]])
    """
    
    # pre-define gradient:
    grad = np.zeros((rays, bins))
    
    # loop over rays:
    for i in range(0, rays):
    
        # only continue for rays with clean bins:
        cleanBins = mask[i][(mask[i]==1)]
        if len(cleanBins) > 0:
            # calculate gradient along ray in (degrees per km):
            grad[i, :] = np.gradient(phidp[i, :], (binLength/1000.0))
        
    # kdp is half of range derivative:
    kdp = 0.5 * grad
    
    # get rid of negative values:
    kdp[(kdp < 0.0)] = 0.0
    
    
    return kdp

# -----------------------------------------------------------------------------
# calculate kdp (method 3 - faster):
# -----------------------------------------------------------------------------

def calc_kdp_v3(rays, bins, binLength, phidp, mask):
    """Calculates KDP as the gradient of reconstructed phidp along each ray
    
    Uses numpy.gradient function - quicker than linear regression

    >>> import numpy as np; calc_kdp_v2(2,2,.6,np.array([[0,6],[0,6]]), np.ones((2,2)))
    info: ignoring rays with no clean pixels ...
    array([[ 5000.,  5000.],
           [ 5000.,  5000.]])
    """
    
    # pre-define gradient:
    grad = np.zeros((rays, bins))
    
    grad[:] = np.gradient(phidp[:], (binLength/1000.0))[1]
        
    # kdp is half of range derivative:
    kdp = 0.5 * grad
    
    # get rid of negative values:
    #kdp[(kdp < 0.0)] = 0.0
    
    return kdp
# -----------------------------------------------------------------------------
# calculate moving average:
# -----------------------------------------------------------------------------

def moving_average(vector, winSize):
    """Along-ray moving average of 1D *vector* with window size *winSize*."""
           
    # smooth using convolve function:
    vector_smooth = np.convolve(vector, np.ones(winSize)/winSize, 'same')
        
    # correct for boundary effects:
    vector_smooth[0:(winSize//2)] = vector_smooth[(winSize//2)]
    vector_smooth[-(winSize//2): ] = vector_smooth[-(winSize//2)-1]
        
    return vector_smooth
    
# -----------------------------------------------------------------------------
# smooth data:
# -----------------------------------------------------------------------------
    
def smooth_data(rays, bins, data, winSize):
    """Smooth 2D array using moving average function for each row."""
    
    # pre-define:
    data_smooth = np.zeros((rays, bins))

    # loop over rays:
    for i in range(0, rays):
        # smooth using moving average function:
        data_smooth[i, :] = moving_average(data[i, :], winSize)
    
    return data_smooth

# -----------------------------------------------------------------------------
# write data:
# -----------------------------------------------------------------------------

def write_data(outFile, pv_apc, kdp):
    """Add new *kdp* field to polar volume *pv_apc*, write all data to 
       *outFile*.
    """
    
    # append kdp field to polar volume: (numpy default type is double, 
    # force to float)
    pv_apc.addData('KDP', kdp.astype(np.float32)) 
     
    # write data to file:
    writer = rn.DataGroupWriteHDF5()
    writer.odim_compliant = False
    writer.Write(outFile, pv_apc) 
                                
    return   
    
# -----------------------------------------------------------------------------
#
# -----------------------------------------------------------------------------

    
    
    















