from __future__ import division
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import pylab as P

def get_observation_significances(file, x_col, y_col):
    '''take the significances file and return a dictionary of all ids and p_values for
    each snps.'''
    result = {}
    y_values = []
    for line in file:
        if line.startswith('OTU'):
            pass
        else:
            fields = line.strip().split('\t')   
            x_value = int(fields[x_col].split(':')[1])
            try:
                y_value = float(fields[y_col])
            except ValueError:
                y_value = 0
            result[x_value] = y_value
    return result
    
def get_bins(observation_significances, bin_size, snps):
    '''Return a list of the bins, the list can be based on either snps or bases'''
    if snps:
        #get a list of keys and values from the dictionary of observations
        bins = observation_significances.items()
        bins.sort()
        #sort this list based on snp location
        bins = [id for id, value in bins]
        #filter the list doen to the bin size. This is based on snps. a bin of 1000 will 
        #contain 1000 snps but the number of total bases will change.
        result = bins[0::bin_size]
        result.append(max(observation_significances))
    else:
        #if the bin size should be based on number of snps the following would create that
        result = range(min(observation_significances),
                 max(observation_significances), 
                 bin_size)
        result.append(max(observation_significances))
    return result

def get_observation_counts(observation_significances, bins, alpha):
    ''' take the list of bins and the dictionary of all snps and return dictionaries. one 
    where the key is each bin and the value is the total number of snps in that bin, and 
    one where the key is the bin and the value is the number of snps with p-values less 
    than alpha.'''
    total_observations = {}
    positive_observations = {}
    bins = np.array(bins)
    #this is import to make sure that every item in bins is present in the dictionaries
    for bin in bins:
        total_observations[bin] = 0
        positive_observations[bin] = 0 
    for id, p_value in observation_significances.items():
        bin = np.max(bins[bins <= id])
        total_observations[bin] += 1
        if p_value <= alpha: 
            positive_observations[bin] += 1
    return total_observations, positive_observations

def get_observation_ratios(file, bin_size, x_col=0, y_col=1, alpha=.05, snps=True):
    '''Return a dictionary where the key is the bin, and the value is the ratio of snps
    with a p-value lower than alpha per expected snps based on alpha.'''
    observation_significances = get_observation_significances(file, 
                                                             x_col, 
                                                             y_col)
    bins = get_bins(observation_significances, 
                    bin_size, snps)
    total_observations, positive_observations = \
    get_observation_counts(observation_significances, 
    bins, 
    alpha)
    result = {}
    for bin in bins:
        expected_pos = total_observations[bin]*alpha
        actual_pos = positive_observations[bin]
        try:
            result[bin] = actual_pos/expected_pos
        except ZeroDivisionError:
            result[bin] = 1 
    return result, bins

        
def main():

    file = open('/Users/jc33/Dropbox/caporaso_lab/vcf_to_biom/biom_analysis/out_sig_22.txt', 'U')
    bin_size = 10000
    alpha = .05
    
    bined_data, bins1 = get_observation_ratios(file, bin_size, alpha=alpha)
    P.hist(bined_data.keys(), weights=bined_data.values(), bins=bins1)
    P.plot([1e7, 5.3e7], [1, 1], 'k-', lw=2, color='red')
    P.xlabel('Genome Location')
    P.ylabel('P_value Ratio')
    P.show()

if __name__ == "__main__":
    main()