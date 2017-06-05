# -*- coding: utf-8 -*-
"""
Created on Thu Jan  5 11:17:32 2017

@author: lg1u16
"""

import numpy as np
import pandas as pd
from scipy.ndimage.filters import generic_filter
import math
from scipy import ndimage

# SET UP FUNCTIONS #############################################################

# function to get binary classification from MPDM
def mpd_prop(nRow, nCol, h, p):
    """Create a neutral landscape which has a specific spatial autocorrelation 
    and proportion
    
    Parameters
    ----------
    nRow: int
        Number of rows in the landscape
    nCol: int
        Number of columns in the landscape
    h: float
        Level of spatial autocorrelation - 0 is random, 1 is totally correlated
    p: float
        Proportion of landcover 1 in the landscape
        
    Returns
    -------
    prop_out: np.array
        NLM with required features
    """
    mpd_out = mpd(nRow, nCol, h)
    prop_out = classifyArray(mpd_out, [1 - p, p])
    return prop_out;
    
# function to get binary classification from MPDM - put into df form
def mpd_prop_df(nRow, nCol, p, h):
    """Create a neutral landscape which has a specific spatial autocorrelation 
    and proportion. This outputs as a dataframe for plotting in R.
    
    Parameters
    ----------
    nRow: int
        Number of rows in the landscape
    nCol: int
        Number of columns in the landscape
    h: float
        Level of spatial autocorrelation - 0 is random, 1 is totally correlated
    p: float
        Proportion of landcover 1 in the landscape
        
    Returns
    -------
    out: pd.DataFrame
        NLM with required features
    """
    
    mpd_out = mpd(nRow, nCol, h)
    prop_out = classifyArray(mpd_out, [1 - p, p])
    out = pd.DataFrame(prop_out)
    out['y'] = range(0, nRow, 1)
    out = pd.melt(out, id_vars=['y'], var_name='x', value_name='value')
    out['h'] = h
    out['p'] = p
    return out;


def dd_simple_esmod(params, nRow, nCol):
    """Simplest ES model where the ES value a distance weighting on the value of 
    the land cover within a particular window. ES1 is currently comparable to 
    the boring one and ES2 is comparable to ES1 in the two step model
    
    Parameters
    ----------
    params: pd.series
        Values of each of the parameters (p, h, sigma - the 'scale' parameter 
        for the distance decay function, r - replicate)
    nRow: int
        Number of rows in the landscape
    nCol: int
        Number of columns in the landscape
    
    """
    p = params['p']
    h = params['h']
    r = params['r']
    sigma = params['sigma']

    # simulated landscape
    out = mpd_prop(nRow, nCol, h, p)
    
    # create the weighting surface which will be included as an input to the function
    x, y = np.meshgrid(np.arange(nRow), np.arange(nRow))
    x, y = (x/(nRow-1), y/(nRow-1))
    centre = (0.5, 0.5)
    #centre = (math.floor(nRow/2), math.floor(nRow/2))
    d = np.sqrt((x-centre[1])**2 + (y-centre[0])**2)
    wt = np.exp(-d**2/2*sigma**2)
    wt = wt/np.sum(wt)
    wt = wt.flatten()
      
    wdw1 = generic_filter(out, dd_func, nRow, mode='wrap', extra_keywords = {'lc':1, 'wt':wt})
    wdw2 = wdw1*out
    # output values
    es1_mean = np.mean(wdw1)
    es1_total = np.sum(wdw1)
    es1_var = np.var(wdw1) # NB this is population variance, try to work out if this is right, if sample variance needed use ddof = 1#
    
    es2_mean = np.mean(wdw2)
    es2_total = np.sum(wdw2)
    es2_var = np.var(wdw2)
    return pd.Series({'p_val': p, 'h_val': h, 'rep': r, 'sigma': sigma, 'es1_mean': es1_mean, 'es1_total': es1_total, 'es1_var': es1_var, 'es2_mean': es2_mean, 'es2_total': es2_total, 'es2_var': es2_var});

# function to bring together creating the landscape and moving window analysis in a two step approach (e.g. pollination as input to agri)
def two_step_esmod(params, nRow, nCol):
    """Two stage ES model where first an intermediary layer is calculated as the 
    same as the simple ES model but only if the land cover type is that required. 
    The ES value is then calculated as the mean of the intermediary layer within 
    the desired window.
    
    Parameters
    ----------
    params: pd.series
        Values of each of the parameters (p, h, w - window size, r - replicate)
    nRow: int
        Number of rows in the landscape
    nCol: int
        Number of columns in the landscape    
    """
    p = params['p']
    h = params['h']
    w = int(params['w'])
    r = params['r']
    out = mpd_prop(nRow, nCol, h, p)
    
    # define the ES function
    # output the first layer (e.g. the one where focal patch must be = 1)
    wdw1 = generic_filter(out, np.mean, w, mode='wrap')
    # multiply wdw by out to set the zeros to zero
    wdw1 = wdw1 * out
    # this is currently set to take the relationship as 1:1 (i.e. 10% natural cover = 10% ecosystem service)
    # will need to add in the relationship to create ES surface at a later date
    es1_mean = np.mean(wdw1)
    es1_total = np.sum(wdw1)
    es1_var = np.var(wdw1) # NB this is population variance, try to work out if this is right, if sample variance needed use ddof = 1
    
    # output the second layer (e.g. the one which takes in the first as input and works on the opposite land cover = 0)
    wdw2 = generic_filter(wdw1, np.mean, w, mode='wrap')
    wdw2 = wdw2 * (1 - out)
    es2_mean = np.mean(wdw2)
    es2_total = np.sum(wdw2)
    es2_var = np.var(wdw2) # NB this is population variance, try to work out if this is right, if sample variance needed use ddof = 1
    
    return pd.Series({'p_val': p, 'h_val': h, 'rep': r, 'window_size': w, 'es1_mean': es1_mean, 'es1_total': es1_total, 'es1_var': es1_var, 'es2_mean': es2_mean, 'es2_total': es2_total, 'es2_var': es2_var});

# create functions for exploring the impact of the shape of the relationship between LC proportion and ES output
def lc_prop(values, lc):
    """Calculate the total proportion of a specified list of land covers.
    
    Parameters
    ----------
    values: array
        The array to process
    lc: int, or list of ints
        The land cover number(s) to calculate the proportion of
    
    Returns
    -------
    out: float
        Proportion of all specified land covers within the values {0, 1}  
    """
    out = np.in1d(values,lc).sum(dtype='float') / values.size
    return out

def shannon(values, lc):
    """Calculate the shannon evefor a specified list of land covers.
    
    Parameters
    ----------
    values: array
        The array to process
    lc: int, or list of ints
        The land cover number(s) to calculate the proportion of
    
    Returns
    -------
    out: float
        Shannon evenness of the specified land covers {0, 1}  
    """
    shannon = 0
    
    if np.in1d(values, lc).sum(dtype='float') == 0:
        H = 0
    else:
        for i in lc:
            p = np.in1d(values,i).sum(dtype='float') / np.in1d(values, lc).sum(dtype='float')
            if p == 0:
                shannon = shannon + 0
            else:
                shannon = shannon + (p * np.log(p))
        
        H = -shannon/np.log(len(lc))
    return H

def exp_func(values, lc, a):
    """Calculate the ES value based on an exponential relationship between the 
    mean LC value and the ES value. The returned value is scaled between 0 and 1.
    
    Parameters
    ----------
    values: array
        The array to process
    lc: int
        The land cover number to calculate the proportion of
    a: float
        The slope of the exponential function - the larger a is, the steeper the slope 
        a > 0 has exponential rise, a < 0 has exponential decline
        
    Returns
    -------
    es_value: float
        This is the value of the ES based on the defined rules  {0, 1}
    """
    lc_prop = int((values == lc).sum()) / values.size
    exp_value = np.exp(a*lc_prop)
    es_value = (exp_value - np.exp(a*0)) / (np.exp(a*1) - np.exp(a*0))
    return es_value

def dd_func(values, lc, wt):
    """Calculate the ES value where the desired LC in each cell is weighted by 
    it's distance to the central cell. 
    
    Parameters
    ----------x
    values: array
        The array to process
    lc: int
        The land cover number to calculate the proportion of
    w: array
        A distance weighted surface
    
    Returns
    -------
    es_value: float
        This is the value of the ES based on the defined rules  {0, 1}
    """
    lc_one = (values == lc)*1
    es_value = (lc_one*wt).sum()
    return es_value
       
def simple_model(values, w):
    lin = generic_filter(values, np.mean, 3, mode='wrap')
    exp = (np.exp(lin*5) - np.exp(5*0)) / (np.exp(5*1) - np.exp(5*0))
    negexp = (np.exp(lin*-5) - np.exp(-5*0)) / (np.exp(-5*1) - np.exp(-5*0))
    shannon = np.nan_to_num(-((lin*np.log(lin))+((1-lin)*np.log((1-lin))))/np.log(2))
    return {'window': w, 'linear': np.mean(lin), 'exp': np.mean(exp), 'negexp': np.mean(negexp), 'shannon': np.mean(shannon)}

def apply_function(dat, fn, window_no, window_size):
    """Apply a given function to each cell of a provided np.array. 
    Any new functions to be investigated should be tested here. 
    
    Parameters
    ----------
    dat: array
        The array to process
    fn: string
        The name of the function to apply
    window_no: int
        Is this the first or second window (this affects the standardisation)
        
    Returns
    -------
    out: array
        Processed array where each cell has had the required function applied to it.
    """
    if(window_no == 1):
        max_val = 1
    if(window_no == 2):
        max_val = ((window_size**2 - 1)/(window_size**2))
    
    if(fn == "linear"):
        out = dat / max_val
    elif(fn == "exp"):
        out = (np.exp(dat*5) - np.exp(5*0)) / (np.exp(5*max_val) - np.exp(5*0))
    elif(fn == "negexp"):
        out = (np.exp(dat*-5) - np.exp(-5*0)) / (np.exp(-5*max_val) - np.exp(-5*0))
    elif(fn == "shannon"):
        out = np.nan_to_num(-((dat*np.log(dat))+((1-dat)*np.log((1-dat))))/np.log(2))
        
    return out

def two_step_binary(ls_size, p, h, w1, w2, f1, f2, fp1_same = True, fp2_same = True):
    # create the landscape based on the input parameters
    ls = mpd_prop(ls_size, ls_size, h, p)
    
    # get the mean value of the 'desired' habitat type
    w1_out = generic_filter(ls, np.mean, w1, mode='wrap')
    if(fp1_same == True):
        w1_out = w1_out * ls # this means focal patch type matters
    
    # apply the correct function
    w1_out = apply_function(w1_out, f1, 1, w1)
    
    # get the mean value of the previous output
    w2_out = generic_filter(w1_out, np.mean, w2, mode='wrap')
    if(fp2_same == True):
        w2_out = w2_out * (1 - ls)
    
    w2_out = apply_function(w2_out, f2, 2, w2)
    
    # create output
    return pd.DataFrame({'ls_size': ls_size, 'p_val': p, 'h_val': h, 'w1': w1, 'w2': w2, 'f1': f1, 'f2': f2, 'fp1_same': fp1_same, 'fp2_same': fp2_same, 'es_mean': np.mean(w2_out), 'es_var': np.var(w2_out)}, index=[0])
    
def one_step_binary(ls_size, p, h, w1, f1, fp1_same = True):
    # create the landscape based on the input parameters
    ls = mpd_prop(ls_size, ls_size, h, p)
    
    # get the mean value of the 'desired' habitat type
    w1_out = generic_filter(ls, np.mean, w1, mode='wrap')
    if(fp1_same == True):
        w1_out = w1_out * ls # this means focal patch type matters
    
    # apply the correct function
    w1_out = apply_function(w1_out, f1)
    
    # create output
    return pd.Series({'ls_size': ls_size, 'p_val': p, 'h_val': h, 'w1': w1, 'f1': f1, 'fp1_same': fp1_same, 'es_mean': np.mean(w1_out), 'es_var': np.var(w1_out)})
  
def two_step_binary_cont(ls_size, p, h1, h2, w1, w2, w3, fp1_same = True, fp2_same = True):
    # note w1 is for the first stage, w2 is the scale at which the output of first
    # stage is important, w3 is for the continuous ls
    
    # at the moment this function doesn't allow for differing functions on the links
    ls_binary = mpd_prop(ls_size, ls_size, h1, p)
    ls_cont = mpd(ls_size, ls_size, h2)
    
    w1_out = generic_filter(ls_binary, np.mean, w1, mode='wrap')
    if(fp1_same == True):
        w1_out = w1_out * ls_binary
        
    w2_out = generic_filter(w1_out, np.mean, w2, mode='wrap')
    if(fp2_same == True):
        w2_out = w2_out * (1 - ls_binary)
        
    w3_out = generic_filter(ls_cont, np.var, w3, mode='wrap') # this is the context dependency 
    
    out = w2_out * (1 - w3_out) # in this instance, the more variable the continuous surface within the window, the less the effect of the pollinators
    
    return pd.Series({'ls_size': ls_size, 'p_val': p, 'h_val1': h1, 'h_val2': h2, 'w1': w1, 'fp1_same': fp1_same, 'w2': w2, 'fp2_same': fp2_same, 'w3': w3, 'es_mean': np.mean(out), 'es_var': np.var(out)})

def farmland_birds_sim(ls_size, h, w1, w2, npp):
    """Function to predict species richness of farmland bird indicator species using amount 
    and heterogeneity of habitat at the appropriate scale. Currently the proportions for each 
    landscape are fixed, only the spatial autocorrelation changes
    
    Parameters
    ----------
    ls_size: int
        Side length of the landscape
    h: float
        spatial autocorrelation of the landscape
    w1: int
        size of the amount window
    w2: int
        size of the heterogeneity window
    npp: float
        level of npp for the landscape
        
    Returns
    -------
    out: array
        Predicted species richness of farmland bird indicator species
    """
    
    # create landscape with four land cover types 1-3 are habitat, 0 is not
    ls = mpd(ls_size, ls_size, h)
    ls = classifyArray(ls, [0.25, 0.25, 0.25, 0.25])
    binary = (ls != 0)*1
    # for each cell, calculate habitat amount within the window
    ls_amount = generic_filter(ls, lc_prop, w1, mode='wrap', extra_keywords = {'lc':[1,2,3]})*binary
    # for each cell, calculate the habitat heterogeneity within the window
    ls_hetero = generic_filter(ls, shannon, w2, mode='wrap', extra_keywords = {'lc':[1,2,3]})
    # multiply the amount*hetero*npp
    out = ls_amount * ls_hetero * npp
    return pd.Series({'ls_size': ls_size, 'h_val': h, 'w1': w1, 'w2': w2, 'npp': npp, 'es_mean': np.mean(out), 'es_var': np.var(out)})
    
    
def get_cell_buffer(grid_dat, grid_ref, buffer_size):
    """Function to extract a specified cell from an input grid and then create 
    a square buffer of a fixed size around it. 
    
    Parameters
    ----------
    grid_dat: geopandas GeoDataFrame
        Grid file for the study area
    grid_ref: string
        Grid reference required
    buffer_size: int
        Size of required buffer - will be ~ half the window size for analysis. 
    
    Returns
    -------
    cellb: geopandas GeoDataSeries
        Selected grid cell with required buffer. 
    """
        
    # select a single bng 10 km cell, add a buffer of (w-1)/2, where w = window size, 
    # and clip the LCM raster by this cell
    cell = grid_dat.query('TILE_NAME == "' + grid_ref + '"')

    # this has been checked and creates expected area
    cellb = cell.geometry.apply(lambda g: g.buffer(buffer_size, cap_style=3, join_style=2))
    
    return(cellb)















# FUNCTIONS AT THE BOTTOM ARE ALL TAKEN FROM THE NLMPY PACKAGE (ETHERINGTON ET AL. ) BECAUSE IT'S EASIER THAN INSTALLING ON THE HPC
#------------------------------------------------------------------------------
# REQUIRED FUNCTIONS:
#------------------------------------------------------------------------------

def linearRescale01(array):
    """    
    A rescale in which the values in the array are linearly rescaled to range
    between 0 and 1.

    Parameters
    ----------
    array : array
        2D array of data values.
        
    Returns
    -------
    out : array
        2D array with rescaled values.
    """
    rescaledArray = (array - np.nanmin(array)) / np.nanmax(array - np.nanmin(array))
    return(rescaledArray)

#------------------------------------------------------------------------------

# A function to insert nan cells into an array based on a binary mask array.
def maskArray(array, maskArray):
    """    
    Return the array with nan values inserted where present in the mask array.
    It is assumed that both the arrays have the same dimensions.

    Parameters
    ----------
    array : array
        2D array of data values.
    maskArray : array
        2D array used as a binary mask.
        
    Returns
    -------
    out : array
        2D array with masked values.
    """
    np.place(array, maskArray==0, np.nan)
    return(array)

#------------------------------------------------------------------------------

def randomUniform01(nRow, nCol, mask=None):
    """    
    Create an array with random values ranging 0-1.

    Parameters
    ----------
    nRow : int
        The number of rows in the array.
    nCol : int
        The number of columns in the array.
    mask : array, optional
        2D array used as a binary mask to limit the elements with values.
        
    Returns
    -------
    out : array
        2D float array.
    """
    if mask is None:
        mask = np.ones((nRow, nCol))
    array = np.random.random((nRow, nCol))
    maskedArray = maskArray(array, mask)
    rescaledArray = linearRescale01(maskedArray)
    return(rescaledArray)
    
#------------------------------------------------------------------------------

def nnInterpolate(array, missing):
    """    
    Two-dimensional array nearest-neighbour interpolation in which the elements
    in the positions indicated by the array "missing" are replaced by the
    nearest value from the "array" of data values.

    Parameters
    ----------
    array : array
        2D array of data values.
    missing: boolean array
        Values of True receive interpolated values.
        
    Returns
    -------
    out : array
        2D array with interpolated values.
    """
    # Get row column based index of nearest value
    rcIndex = ndimage.distance_transform_edt(missing, return_distances=False, 
                                             return_indices=True)
    # Create a complete array by extracting values based on the index
    interpolatedArray = array[tuple(rcIndex)]
    return(interpolatedArray)

#------------------------------------------------------------------------------

def w2cp(weights):
    """    
    Convert a list of category weights into a 1D NumPy array of cumulative 
    proportions.

    Parameters
    ----------
    weights : list
        A list of numeric values
        
    Returns
    -------
    out : array
        1D array of class cumulative proportions.
    """
    w = np.array(weights, dtype=float)
    proportions = w / np.sum(w)
    cumulativeProportions = np.cumsum(proportions)
    cumulativeProportions[-1] = 1 # to ensure the last value is 1
    return(cumulativeProportions)

#------------------------------------------------------------------------------

def calcBoundaries(array, cumulativeProportions, classifyMask=None):
    """    
    Determine upper class boundaries for classification of an array with values
    ranging 0-1 based upon an array of cumulative proportions.

    Parameters
    ----------
    array : array
        2D array of data values.
    cumulativeProportions : array
        1D array of class cumulative proportions.
    classifyMask : array, optional
        2D array used as a binary mask to limit the elements used to determine
        the upper boundary values for each class.
        
    Returns
    -------
    out : array
        1D float array.
    """
    if classifyMask is None:
        classifyMask = np.ones(np.shape(array))
    maskedArray = array * classifyMask
    np.place(maskedArray, classifyMask==0, np.nan)
    # Determine the number of cells that are in the classification mask.
    nCells = np.count_nonzero(np.isfinite(maskedArray))
    # Based on the number of cells, find the index of upper boundary element
    boundaryIndexes = (cumulativeProportions * nCells).astype(int) - 1
    # Index out the the upper boundary value for each class
    boundaryValues = np.sort(np.ndarray.flatten(maskedArray))[boundaryIndexes]
    # Ensure the maximum boundary value is equal to 1
    boundaryValues[-1] = 1
    return(boundaryValues)

#------------------------------------------------------------------------------
      
def classifyArray(array, weights, classifyMask=None):
    """    
    Classify an array with values ranging 0-1 into proportions based upon a 
    list of class weights.

    Parameters
    ----------
    array : array
        2D array of data values.
    weights : list
        A list of numeric values
    classifyMask : array, optional
        2D array used as a binary mask to limit the elements used to determine
        the upper boundary values for each class.
        
    Returns
    -------
    out : array
        2D array.
    """
    cumulativeProportions = w2cp(weights)
    boundaryValues = calcBoundaries(array, cumulativeProportions, classifyMask)
    # Classify the array
    classifiedArray = np.searchsorted(boundaryValues, array)
    # Replace any nan values
    classifiedArray = classifiedArray.astype(float)
    np.place(classifiedArray, np.isnan(array), np.nan)
    return(classifiedArray)

#------------------------------------------------------------------------------

def blendArray(primaryArray, arrays, scalingFactors=None):
    """    
    Blend a primary array with other arrays weighted by scaling factors.

    Parameters
    ----------
    primaryArray : array
        2D array of data values.
    arrays : list
        List of 2D arrays of data values.
    scalingFactors : list
        List of scaling factors used to weight the arrays in the blend.
        
    Returns
    -------
    out : array
        2D array.
    """
    if scalingFactors is None:
        scalingFactors = np.ones(len(arrays))
    for n in range(len(arrays)):
        primaryArray = primaryArray + (arrays[n] * scalingFactors[n])
    blendedArray = primaryArray / len(arrays)
    rescaledArray = linearRescale01(blendedArray)
    return(rescaledArray)
    
#------------------------------------------------------------------------------

def blendClusterArray(primaryArray, arrays, scalingFactors=None):
    """    
    Blend a primary cluster NLM with other arrays in which the mean value per 
    cluster is weighted by scaling factors.

    Parameters
    ----------
    primaryArray : array
        2D array of data values in which values are clustered.
    arrays : list
        List of 2D arrays of data values.
    scalingFactors : list
        List of scaling factors used to weight the arrays in the blend.
        
    Returns
    -------
    out : array
        2D array.
    """
    if scalingFactors is None:
        scalingFactors = np.ones(len(arrays))
    for n in range(len(arrays)):
        meanOfClusterArray = meanOfCluster(primaryArray, arrays[n])
        primaryArray = primaryArray + (meanOfClusterArray * scalingFactors[n])
    blendedArray = primaryArray / len(arrays)
    rescaledArray = linearRescale01(blendedArray)
    return(rescaledArray)
    
#------------------------------------------------------------------------------

def meanOfCluster(clusterArray, array):
    """    
    For each cluster of elements in an array, calculate the mean value for the
    cluster based on a second array.

    Parameters
    ----------
    clutserArray : array
        2D array of data values in which values are clustered.
    array : array
        2D array of data values.
        
    Returns
    -------
    out : array
        2D array.
    """
    meanClusterValues = np.zeros(np.shape(clusterArray))
    clusterValues = np.unique(clusterArray)
    for value in clusterValues:
        if np.isfinite(value):
            # Extract location of values
            valueLocs = clusterArray == value
            # Define clusters in array
            clusters, nClusters = ndimage.measurements.label(valueLocs)
            # Get mean for each cluster
            means = ndimage.mean(array, clusters, range(1,nClusters + 1))
            means = np.insert(means, 0, 0) # for background non-cluster
            # Apply mean values to clusters by index
            clusterMeans = means[clusters]
            # Add values for those clusters to array
            meanClusterValues = meanClusterValues + clusterMeans
    np.place(meanClusterValues, np.isnan(clusterArray), np.nan)
    rescaledArray = linearRescale01(meanClusterValues)
    return(rescaledArray)

#------------------------------------------------------------------------------

def exportASCIIGrid(outFile, nlm, xll=0, yll=0, cellSize=1):
    """
    Export a NLM array as a ASCII grid raster file.
    
    Parameters
    ----------
    outFile : string
        The path and name of the output raster file.
    nlm : 2D array
        The NLM to be exported.
    xll : number
        Raster lower left corner x coordinate.
    yll : number
        Raster lower left corner y coordinate.
    cellSize : number
        The size of the cells in the output raster.
    """
    # Get dimensions of the NLM
    nRow, nCol = nlm.shape
    # Convert any nan elements to null data value of -9999
    np.place(nlm, np.isnan(nlm), -9999)
    # Create raster out file
    textOut = open(outFile, 'w')
    # Write metadata
    textOut.write("NCOLS " + str(nCol) + "\n")
    textOut.write("NROWS " + str(nRow) + "\n")
    textOut.write("XLLCORNER " + str(xll) + "\n")
    textOut.write("YLLCORNER " + str(yll) + "\n")
    textOut.write("CELLSIZE " + str(cellSize) + "\n")
    textOut.write("NODATA_VALUE -9999\n")
    # Write NLM
    for row in range(nRow):
        lineout = ""
        for col in range(nCol):
            lineout = lineout + str(nlm[row,col]) + " "
        textOut.write(lineout[:-1] + "\n")
    textOut.close()

#------------------------------------------------------------------------------
# NEUTRAL LANDSCAPE MODELS:
#------------------------------------------------------------------------------

def random(nRow, nCol, mask=None):
    """    
    Create a spatially random neutral landscape model with values ranging 0-1.

    Parameters
    ----------
    nRow : int
        The number of rows in the array.
    nCol : int
        The number of columns in the array.
    mask : array, optional
        2D array used as a binary mask to limit the elements with values.
        
    Returns
    -------
    out : array
        2D array.
    """
    array = randomUniform01(nRow, nCol, mask)
    return(array)
    
#------------------------------------------------------------------------------

def planarGradient(nRow, nCol, direction=None, mask=None):
    """    
    Create a planar gradient neutral landscape model with values ranging 0-1.

    Parameters
    ----------
    nRow : int
        The number of rows in the array.
    nCol : int
        The number of columns in the array.
    direction: int, optional
        The direction of the gradient as a bearing from north, if unspecified
        the direction is randomly determined.
    mask : array, optional
        2D array used as a binary mask to limit the elements with values.
        
    Returns
    -------
    out : array
        2D array.
    """
    if direction is None:
        direction = np.random.uniform(0, 360, 1) # a random direction
    if mask is None:
        mask = np.ones((nRow, nCol))
    # Create arrays of row and column index
    rowIndex, colIndex = np.indices((nRow, nCol))
    # Determine the eastness and southness of the direction
    eastness = np.sin(np.deg2rad(direction))
    southness = np.cos(np.deg2rad(direction)) * -1
    # Create gradient array
    gradientArray = (southness * rowIndex + eastness * colIndex)
    maskedArray = maskArray(gradientArray, mask)
    rescaledArray = linearRescale01(maskedArray)
    return(rescaledArray)

#------------------------------------------------------------------------------

def edgeGradient(nRow, nCol, direction=None, mask=None):
    """    
    Create an edge gradient neutral landscape model with values ranging 0-1.

    Parameters
    ----------
    nRow : int
        The number of rows in the array.
    nCol : int
        The number of columns in the array.
    direction: int, optional
        The direction of the gradient as a bearing from north, if unspecified
        the direction is randomly determined.
    mask : array, optional
        2D array used as a binary mask to limit the elements with values.
        
    Returns
    -------
    out : array
        2D array.
    """
    # Create planar gradient
    gradientArray = planarGradient(nRow, nCol, direction, mask)
    # Transform to a central gradient
    edgeGradientArray = (np.abs(0.5 - gradientArray) * -2) + 1
    rescaledArray = linearRescale01(edgeGradientArray)
    return(rescaledArray)

#------------------------------------------------------------------------------

def distanceGradient(source, mask=None):
    """    
    Create a distance gradient neutral landscape model with values ranging 0-1.

    Parameters
    ----------
    source : array
        2D array binary array that defines the source elements from which
        distance will be measured.  The dimensions of source also specify
        the output dimensions of the distance gradient.
    mask : array, optional
        2D array used as a binary mask to limit the elements with values.
        
    Returns
    -------
    out : array
        2D array.
    """
    if mask is None:
        mask = np.ones(np.shape(source))
    gradient = ndimage.distance_transform_edt(1 - source)
    maskedArray = maskArray(gradient, mask)
    rescaledArray = linearRescale01(maskedArray)
    return(rescaledArray)

#------------------------------------------------------------------------------

def mpd(nRow, nCol, h, mask=None):
    """    
    Create a midpoint displacement neutral landscape model with values ranging 
    0-1.

    Parameters
    ----------
    nRow : int
        The number of rows in the array.
    nCol : int
        The number of columns in the array.
    h: float
        The h value controls the level of spatial autocorrelation in element
        values.
    mask : array, optional
        2D array used as a binary mask to limit the elements with values.
        
    Returns
    -------
    out : array
        2D array.
    """
    if mask is None:
        mask = np.ones((nRow, nCol))  
    # Determine the dimension of the smallest square
    maxDim = max(nRow, nCol)
    N = int(math.ceil(math.log(maxDim - 1, 2)))
    dim = 2 ** N + 1
    # Create a surface consisting of random displacement heights average value
    # 0, range from [-0.5, 0.5] x displacementheight
    disheight = 2.0
    surface = np.random.random([dim,dim]) * disheight -0.5 * disheight
    
    #--------------------------------------------------------------------------
    
    # Apply the square-diamond algorithm
    def randomdisplace(disheight):
        # Returns a random displacement between -0.5 * disheight and 0.5 * disheight
        return np.random.random() * disheight -0.5 * disheight
    
    def displacevals(p, disheight):
        # Calculate the average value of the 4 corners of a square (3 if up
        # against a corner) and displace at random.
        if len(p) == 4:
            pcentre = 0.25 * sum(p) + randomdisplace(disheight)
        elif len(p) == 3:
            pcentre = sum(p) / 3 + randomdisplace(disheight)	
        return pcentre
    
    def check_diamond_coords(diax,diay,dim,i2):
        # get the coordinates of the diamond centred on diax, diay with radius i2
        # if it fits inside the study area
        if diax < 0 or diax > dim or diay <0 or diay > dim:
            return []
        if diax-i2 < 0:
            return [(diax+i2,diay),(diax,diay-i2),(diax,diay+i2)]
        if diax + i2 >= dim:
            return [(diax-i2,diay),(diax,diay-i2),(diax,diay+i2)]
        if diay-i2 < 0:
            return [(diax+i2,diay),(diax-i2,diay),(diax,diay+i2)]
        if diay+i2 >= dim:
            return [(diax+i2,diay),(diax-i2,diay),(diax,diay-i2)]
        return [(diax+i2,diay),(diax-i2,diay),(diax,diay-i2),(diax,diay+i2)]

    # Set square size to cover the whole array
    inc = dim-1
    while inc > 1: # while considering a square/diamond at least 2x2 in size
            
            i2 = int(inc/2) # what is half the width (i.e. where is the centre?)
            # SQUARE step
            for x in range(0,dim-1,inc):
                    for y in range(0,dim-1,inc):
                            # this adjusts the centre of the square 
                            surface[x+i2,y+i2] = displacevals([surface[x,y],surface[x+inc,y],surface[x+inc,y+inc],surface[x,y+inc]],disheight)
            
            # DIAMOND step
            for x in range(0, dim-1, inc):
                for y in range(0, dim-1,inc):
                    diaco = check_diamond_coords(x+i2,y,dim,i2)
                    diavals = []
                    for co in diaco:
                        diavals.append(surface[co])
                    surface[x+i2,y] = displacevals(diavals,disheight)
                   
                    diaco = check_diamond_coords(x,y+i2,dim,i2)
                    diavals = []
                    for co in diaco:
                        diavals.append(surface[co])
                    surface[x,y+i2] = displacevals(diavals,disheight)

                    diaco = check_diamond_coords(x+inc,y+i2,dim,i2)
                    diavals = []
                    for co in diaco:
                        diavals.append(surface[co])
                    surface[x+inc,y+i2] = displacevals(diavals,disheight)

                    diaco = check_diamond_coords(x+i2,y+inc,dim,i2)
                    diavals = []
                    for co in diaco:
                        diavals.append(surface[co])
                    surface[x+i2,y+inc] = displacevals(diavals,disheight)
                    
            # Reduce displacement height
            disheight = disheight * 2 ** (-h)
            inc = int(inc / 2)

    #--------------------------------------------------------------------------
    
    # Extract a portion of the array to match the dimensions
    randomStartRow = np.random.choice(range(dim - nRow))
    randomStartCol = np.random.choice(range(dim - nCol))
    array = surface[randomStartRow:randomStartRow + nRow,
                    randomStartCol:randomStartCol + nCol]
    # Apply mask and rescale 0-1
    maskedArray = maskArray(array, mask)
    rescaledArray = linearRescale01(maskedArray)
    return(rescaledArray)
    
#------------------------------------------------------------------------------

def randomRectangularCluster(nRow, nCol, minL, maxL, mask=None):
    """    
    Create a random rectangular cluster neutral landscape model with 
    values ranging 0-1.

    Parameters
    ----------
    nRow : int
        The number of rows in the array.
    nCol : int
        The number of columns in the array.
    minL: int
        The minimum possible length of width and height for each random 
        rectangular cluster.
    maxL: int
        The maximum possible length of width and height for each random 
        rectangular cluster.
    mask : array, optional
        2D array used as a binary mask to limit the elements with values.
        
    Returns
    -------
    out : array
        2D array.
    """    
    if mask is None:
        mask = np.ones((nRow, nCol))
    # Create an empty array of correct dimensions
    array = np.zeros((nRow, nCol)) - 1
    # Keep applying random clusters until all elements have a value
    while np.min(array) == -1:
        width = np.random.choice(range(minL, maxL))
        height = np.random.choice(range(minL, maxL))
        row = np.random.choice(range(-maxL, nRow))
        col = np.random.choice(range(-maxL, nCol))
        array[row:row + width, col:col + height] = np.random.random()   
    # Apply mask and rescale 0-1        
    maskedArray = maskArray(array, mask)
    rescaledArray = linearRescale01(maskedArray)
    return(rescaledArray)

#------------------------------------------------------------------------------

def randomElementNN(nRow, nCol, n, mask=None):
    """    
    Create a random element nearest-neighbour neutral landscape model with 
    values ranging 0-1.

    Parameters
    ----------
    nRow : int
        The number of rows in the array.
    nCol : int
        The number of columns in the array.
    n: int
        The number of elements randomly selected to form the basis of
        nearest-neighbour clusters.
    mask : array, optional
        2D array used as a binary mask to limit the elements with values.
        
    Returns
    -------
    out : array
        2D array.
    """
    if mask is None:
        mask = np.ones((nRow, nCol))
    # Create an empty array of correct dimensions
    array = np.zeros((nRow, nCol))
    # Insert value for n elements
    for element in range(n):
        randomRow = np.random.choice(range(nRow))
        randomCol = np.random.choice(range(nCol))
        if array[randomRow, randomCol] == 0 and mask[randomRow, randomCol] == 1:
            array[randomRow, randomCol] = np.random.random(1)
    # Interpolate the values
    interpolatedArray = nnInterpolate(array, array==0)
    # Apply mask and rescale 0-1
    maskedArray = maskArray(interpolatedArray, mask)
    rescaledArray = linearRescale01(maskedArray)
    return(rescaledArray)
 
#------------------------------------------------------------------------------

def randomClusterNN(nRow, nCol, p, n='4-neighbourhood', mask=None):
    """    
    Create a random cluster nearest-neighbour neutral landscape model with 
    values ranging 0-1.

    Parameters
    ----------
    nRow : int
        The number of rows in the array.
    nCol : int
        The number of columns in the array.
    p: float
        The p value controls the proportion of elements randomly selected to
        form clusters.
    n: string, optional
        Clusters are defined using a set of neighbourhood structures that 
        include:
                            [0,1,0]
        '4-neighbourhood' = [1,1,1]
                            [0,1,0]
                            
                            [1,1,1]
        '8-neighbourhood' = [1,1,1]
                            [1,1,1]
                            
                     [0,1,1]
        'diagonal' = [1,1,1]
                     [1,1,0]
                     
        The default parameter setting is '4-neighbourhood'.
        
    mask : array, optional
        2D array used as a binary mask to limit the elements with values.
        
    Returns
    -------
    out : array
        2D array.
    """
    if mask is None:
        mask = np.ones((nRow, nCol))
    # Define a dictionary of possible neighbourhood structures:
    neighbourhoods = {}
    neighbourhoods['4-neighbourhood'] = np.array([[0,1,0],
                                                  [1,1,1],
                                                  [0,1,0]])
    neighbourhoods['8-neighbourhood'] = np.array([[1,1,1],
                                                  [1,1,1],
                                                  [1,1,1]])
    neighbourhoods['diagonal'] = np.array([[0,1,1],
                                           [1,1,1],
                                           [1,1,0]])
    # Create percolation array
    randomArray = random(nRow, nCol, mask)
    percolationArray = classifyArray(randomArray, [1 - p, p])
    # As nan not supported in cluster algorithm replace with zeros
    np.place(percolationArray, np.isnan(percolationArray), 0)
    # Define clusters
    clusters, nClusters = ndimage.measurements.label(percolationArray, 
                                                     neighbourhoods[n])
    # Create random set of values for each the clusters
    randomValues = np.random.random(nClusters)
    randomValues = np.insert(randomValues, 0, 0) # for background non-cluster
    # Apply values by indexing by cluster
    clusterArray = randomValues[clusters]
    # Gap fill with nearest neighbour interpolation
    interpolatedArray = nnInterpolate(clusterArray, clusterArray==0)
    # Apply mask and rescale 0-1
    maskedArray = maskArray(interpolatedArray, mask)
    rescaledArray = linearRescale01(maskedArray)
    return(rescaledArray)

#------------------------------------------------------------------------------
