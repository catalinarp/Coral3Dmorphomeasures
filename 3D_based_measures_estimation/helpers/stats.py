import numpy as np
import scipy.stats as spy

def getDistributionCharactertics(raw_data, perc_or_z = 'perc', perc=2.5, Z = 3,\
        weighted = False, weights = None, selection=[]):

    # make a selection of the data
    if len(selection) > 0:
        raw_data = raw_data[selection]

    # remove outliers by eihter percentile or Z-score
    if perc_or_z == 'perc':
        outliers = outliersByPerc(raw_data,perc=perc)
    elif perc_or_z == 'Z':
        outliers = outliersByZ(raw_data, Z_thres=Z)
    else:
        print('perc_or_z should be either "Z" or "perc"')

    results = np.delete(raw_data, outliers)

    if weighted:
        wts = np.delete(weights, outliers)
        wts = wts/np.sum(wts)
        
        mean = weightedMean(results,wts)
        var = weightedVariance(results, wts, mean)
        skew = weightedSkew(results, wts, mean, var)
        kurt = weightedKurtosis(results, wts, mean, var)
    else:
        mean = np.mean(results)
        var = np.var(results)
        skew = spy.skew(results)
        kurt = spy.kurtosis(results)

    return {'mean': mean, 'variance': var, 'skewness': skew, 'kurtosis': kurt}

def getZscore(val):
    """ Calculated Z-score """
    Z = (val - np.mean(val))/np.std(val)
    return Z

def weightedMean(val, wts):
    """Calculates the weighted mean"""
    return np.average(val, weights=wts)


def weightedVariance(val, wts, wmean = None):
    """Calculates the weighted variance"""
    
    if not wmean:
        wmean = weightedMean(val,wts)
    return np.average((val - wmean)**2, weights=wts)


def weightedSkew(val, wts, wmean = None, wvar = None):
    """Calculates the weighted skewness"""
    
    if not wmean:
        wmean = weightedMean(val, wrts)
    if not wvar:
        wvar = weightedVariance(val, wrts, wmean)
    
    return (np.average((val - wmean)**3, weights=wts) /
            (wvar**1.5))

def weightedKurtosis(val, wts,wmean = None, wvar = None):
    """Calculates the weighted skewness"""
    if not wmean:
        wmean = weightedMean(val, wrts)
    if not wvar:
        wvar = weightedVariance(val, wrts, wmean)
    
    return (np.average((val - wmean)**4, weights=wts) /
            (wvar**2))

def outliersByZ(val, Z_thres = 3):
    """Identifies outliers based on Z-score"""
    
    Z = getZscore(val)
    Z_abs = np.abs(Z)
    return np.where(Z_abs > Z_thres)

def outliersByPerc(val, perc = 2.5):
    """Identifies outliers based on percentile"""
    
    min_val = np.percentile(val, perc)
    max_val = np.percentile(val, 100 - perc)
    return np.where((val < min_val) | (val > max_val))