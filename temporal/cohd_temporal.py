from copy import deepcopy
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

class AgeCounter():
    def __init__(self, max_age=90):
        self.max_age = max_age
        self.counts = np.zeros((max_age + 1), dtype=np.uint32)
        
    def add_age(self, age):
        self.counts[min(int(age), self.max_age)] += 1
        
    def x(self):
        return np.arange(self.max_age + 1)
        
class GroupedAgeCounts():
    """ max_age: largest age allowed to start a bin """
    def __init__(self, age_counter, bin_width=5, max_age=90):
        self.bin_width = bin_width
        self.max_age = max_age
        self.bins = int(np.ceil((self.max_age + 1) / self.bin_width))
        self.counts = np.zeros(self.bins, dtype=np.uint32)
        
        for i in range(self.bins - 1):
            self.counts[i] = np.sum(age_counter.counts[(i*self.bin_width):((i+1)*self.bin_width)])
        self.counts[self.bins-1] = np.sum(age_counter.counts[((self.bins-1)*self.bin_width):])
        
    def x(self, incremental=False):
        if incremental:
            # Return [0, 1, 2, 3, ...] (good for x-values of bar plot)
            return np.arange(self.bins)
        else:
            # Return e.g. [0, 5, 10, 15, ...] indicating the start of each bin
            return np.arange(self.bins) * self.bin_width
        
    def xlabels(self):
        return [str(x) for x in self.x(incremental=True)*self.bin_width]

class DeltaCounter():
    """Class for keeping counts where the bins grow as power of 2
    
    For example, with n=5, the bins will be
    <=-16 | -15 — -8 | -7 — -4 | -3 — -2 | -1 | 0 | 1 | 2 — 3 | 4 — 7 | 8 — 15 | >= 16
    """
    # Maximum supported n
    _max_n = 14
    
    # Lookup array for quickly finding out the index offset into the counts list    
    _offset = [0] + [(int(np.log2(x)) + 1) for x in range(1, 2**(_max_n - 1))]    
    
    def __init__(self, n=13):
        assert n >= 0
        self.n = min(n, DeltaCounter._max_n)
        self.bins = self.n * 2 + 1  # Number of bins
        
        # The delta at the start of the last bin. delta >= max_delta all go into the last bin
        self.max_delta = 2**(self.n - 1)  
        
        # counts will be indexed from -n:1:n
        self.counts = [0] * self.bins
        
    def add(self, delta):
        if delta == 0:
            self.counts[self.n] += 1
            return
        
        # Determine the index into counts from the delta
        if delta < 0:
            if delta <= -self.max_delta:
                index = 0
            else:
                index = self.n - min(DeltaCounter._offset[-delta], self.n)
        else:
            if delta >= self.max_delta:
                index = self.bins - 1
            else:
                index = self.n + min(DeltaCounter._offset[delta], self.n)
                
        self.counts[index] += 1
        
    def bin_labels(self):
        labels = [''] * self.bins
        labels[self.n] = '0'
        
        if self.n > 0:
            labels[0] = f'-{2**(self.n-1)}+'
            labels[self.bins - 1] = f'{2**(self.n-1)}+'
            labels[slice(self.n - 1, self.n + 2, 2)] = ['-1', '1']
        
        for i in range(2, self.n):
            lower = 2**(i-1)
            upper = 2**i - 1
            labels[self.n + i] = f'{lower} — {upper}'
            labels[self.n - i] = f'-{upper} — -{lower}'
            
        return labels
    
    def bin_labels_mixed(self):
        """ Compact labels with mixed units (easier to interpret) """
        labels = ['0d']
        
        for i in range(0, self.n):
            days = 2**i
            if days >= 365:
                years = days / 365.25
                label = f'{years:0.1f}y'
            else:
                label = f'{days}d'
            
            labels.append(label)
            labels.insert(0, '-' + label)
        
        return labels
    
    def x(self):
        """ X-ticks ranging from -self.n:self.n"""
        return np.arange(-self.n, self.n+1)
    
    def normalized_counts(self):
        """ Normalize the counts relative to the number of days in the bin """
        norm = 2**np.maximum(abs(self.x())-1, 0)
        return self.counts / norm
    
    def get(self, ix):
        """ Gets the counts using a relative index, where 0 is the 0-day bin, 1 is the 1-day bin, 
        -1 is the -1-day bin, etc. 
        
        Params
        ------
        indices: relative index or list of indices """
        if isinstance(ix, int):
            return self.counts[self.n + ix]
        else:
            return [self.counts[self.n + i] for i in ix]
        
    
def get_delta(deltas, concept_source, concept_target):
    """ Retrieves the deltas for concept_source -> concept_target 
    Positive deltas indicate that concept_source occurred before concept_target. 
    Likewise, negative deltas indicate that concept_source occurred after concept_target.

    Params
    ------
    concept_source: the source concept in the pair
    concept_target: the target concept in the pair

    Return
    ------
    The deltas from concept_source -> concept_target. 
    None, if the concept pair did not co-occur or if their total co-occurrence count was too low
    """
    key = (concept_source, concept_target)
    reverse = concept_source > concept_target
    if reverse:
        key = (concept_target, concept_source)
        
    if key in deltas:
        delta = deltas[key]
        if reverse:
            # Create a copy of the DeltaCounter so that reversing the counts doesn't affect the original
            delta = deepcopy(delta)
            delta.counts.reverse()  # Reverses the list in-place (doesn't create a new list)
        return delta
    else:
        return None
    
    
def plot_delta(cohd, concept_source, concept_target, mode='count', bin_normalize=False, alpha=1.0, show_plot=True):
    """ Plots the delta counts from concept_source --> concept_target
    
    mode: 'count' (default), 'density', or 'relative'"""
    cads = cohd.cads
    deltas = cohd.deltas
    
    dc = get_delta(deltas, concept_source, concept_target)
    if dc is None:
        # concept pair did not co-occur or was removed for low count
        return
    x = dc.x() 
                    
    if bin_normalize:
        # normalize the counts relative to the size of the bins
        counts = dc.normalized_counts()
    else:
        # counts not normalized to bin size
        counts = dc.counts    
    
    
    mode = mode.lower()
    if mode == 'density':
        # normalize the co-occurrence counts to a sum of 1.0 (i.e., probability distribution)
        counts = counts / np.sum(counts)
    elif mode == 'relative':
        # normalize by the total counts of the base concept
        counts = counts / np.sum(cads[concept_source].counts)        
        
    plt.bar(x, counts, alpha=alpha)
    plt.xticks(x, dc.bin_labels_mixed(), rotation='vertical')
    concept_source_name, concept_target_name = cohd.concepts.loc[[concept_source, concept_target], 'concept_name']
    plt.title(f'{concept_source_name} -> {concept_target_name}')
    
    if show_plot:
        # Optionally not showing the plot now gives caller ability to edit plot properties
        plt.show()
    
    
def plot_age_distributions(cohd, concept_ids, normalize=True, alpha=None, show_plot=True):
    cads = cohd.cads
    df_concepts = cohd.concepts
    
    if isinstance(concept_ids, int):
        concept_ids = [concept_ids]
        
    if not alpha:
        alpha = 1 if len(concept_ids) == 1 else 0.5     
    
    labels = []
    for c in concept_ids:
        cad = cads[c]
        h = cad.counts
        if normalize:
            h = h / np.sum(h)
        plt.bar(np.arange(cad.max_age + 1), h, alpha=alpha)
        labels.append(df_concepts.loc[c, 'concept_name'])
    
    plt.legend(labels)
    
    if show_plot:
        # Optionally not showing the plot now gives caller ability to edit plot properties
        plt.show()
        

def plot_grouped_age_distribution(gcad, show_plot=True):
    plt.bar(gcad.x(incremental=True), gcad.counts)
    plt.xticks(gcad.x(incremental=True), gcad.xlabels())
    if show_plot:
        plt.show()       
        
        
def plot_delta_simple(delta, show_plot=True):
    x = delta.x()
    counts = delta.counts
    plt.bar(x, counts)
    plt.xticks(x, delta.bin_labels_mixed(), rotation='vertical')
    
    if show_plot:
        # Optionally not showing the plot now gives caller ability to edit plot properties
        plt.show()
        
        
def association_strength(cohd, c1, c2, verbose=False):
    """ Calculates the association strength (i.e., ln_ratio) between c1 and c2 
    
    Return
    ------
    If there are delta counts for the concept pair, returns the association strenght. Else, None."""
    delta = get_delta(deltas, c1, c2)
    if delta is None:
        return None
    
    n_pair = np.sum(delta.counts)
    
    n_c1 = np.sum(cohd.cads[c1].counts)
    n_c2 = np.sum(cohd.cads[c2].counts)
    n_est = (n_c1 * n_c2) / n_patients
    ln_ratio = np.log(n_pair / n_est)
    
    if verbose:        
        print(f'n_c1: {n_c1}\tn_c2: {n_c2}\tn_pair: {n_pair}\tn_est: {n_est}\tln_ratio: {ln_ratio}')
        
    return ln_ratio
    
    
def find_similar_concepts(cohd, concept_id, n=None, method='Jaccard', threshold=None, exc_associated=True,
                          restrict_domain=False):    
    """ Return list of concepts with similar age distributions with the age distribution of concept_id 
    
    Params
    ------
    cohd: COHD data structure
    n: maximum number of similar concepts 
    method: 'Jaccard' (Jaccard similarity) or 'KL' (KL-divergence)
    threshold: (float) Threshold to apply for the similarity metric, default: no threshold applied
    exc_associated: True - exclude similar concepts if they have elevated co-occurrence with the concept of interest
    deltas: The deltas data structure, required if exc_associated is True
    n_patients: 
    restrict_domain: True - only include concepts from the same domain (df_concepts will be required)
    """
    cads = cohd.cads
    
    method = method.strip().lower()
    
    # Get the domain of the concept of interest
    restrict_domain = restrict_domain and df_concepts is not None
    if restrict_domain:
        domain = df_concepts.loc[concept_id, 'domain_id']
        
    # Get concept_id's normalized histogram
    h1 = cads[concept_id].counts
    nh1 = h1 / np.sum(h1)
    
    # Look through age distributions of all concepts to find a similar concept age distribution
    similar_concepts = list()
    metrics = list()
    for c, cad in cads.items():
        # Don't compare against self
        if c == concept_id:
            continue

        # Find concepts of the same domain
        if restrict_domain and df_concepts.loc[c, 'domain_id'] != domain:
            continue

        # Get normalized histogram of concept 2
        h2 = cad.counts
        nh2 = h2 / np.sum(h2)
            
        if method == 'jaccard':
            # Calculate Jaccard similarity
            metric = np.sum(np.minimum(nh1, nh2)) / np.sum(np.maximum(nh1, nh2))
            
            if threshold and metric < threshold:
                # This concept doesn't meet threshold, don't add
                continue
                
        elif method == 'kl':
            # Calculate KL divergence
            metric = scipy.stats.entropy(nh1, nh2)

            # Save concepts that are more similar than the threshold
            if threshold and metric > threshold:
                # This concept doesn't meet threshold, don't add
                continue
                
            similar_concepts.append(c) 
            metrics.append(metric)
            
        # Exclude similar concepts if they are associated with the concept of interest
        if exc_associated:
            # Exclude the similar concept if it has a high association strength with the concept
            ln_ratio = association_strength(cohd, concept_id, c)
            if ln_ratio is None or ln_ratio >= 2.0:                
                continue
                                
            # Exclude the concept if it has relatively high 0-day co-occurrence
            delta = get_delta(deltas, concept_id, c)
            pair_count = np.sum(delta.counts)
            if pair_count < 100:
                # If the co-occurrence count is low, exclude similar concepts if the 0-day count is 1+
                if delta.get(0) > 1:
                    continue
            else:               
                # If the co-occurrence count is high, exclude similar concepts if the 0-day count is more
                # than 1% of the total co-occurrence count
                if (delta.get(0) / pair_count) > 0.01:
                    continue            
            
        similar_concepts.append(c)
        metrics.append(metric)
            
    # Sort the concepts in descending order of similarity
    argsort = np.argsort(metrics)
    if method == 'jaccard':
        argsort = np.flipud(argsort)  # reverse the order
    metrics = [metrics[x] for x in argsort]
    similar_concepts = [similar_concepts[x] for x in argsort]
    
    if n:
        # Return only the top n most similar concepts
        n = min(n, len(metrics))
        metrics = metrics[:n]
        similar_concepts = similar_concepts[:n]
            
    return similar_concepts, metrics
        

def AdaptiveGroupedAgeCountOld(cad, bin_widths=[1, 2, 4, 8, 16], max_age=90, low_count=10, 
                            upper_percentile=90, small_percent=0.04, verbose=False):
    """ Adaptively adjust the grouped concept-age distribution bin width 
    
    The idea is to increase the bin widths until the number of bins with low_count (<=10) is smaller than
    the number of bins that have a small proportion of the counts (less than small_percent% of the upper_percentile)"""
    for bw in bin_widths:
        gcad = GroupedAgeCounts(cad, bin_width=bw, max_age=max_age)
        max_count = np.percentile(gcad.counts, upper_percentile)
        n_less_10 = np.sum(gcad.counts <= low_count)
        n_small_percent = np.sum(gcad.counts <= max_count * small_percent)
        
        if verbose:
            print(f'bin_width: {bw}\tn_less_10: {n_less_10}\tn_small_percent: {n_small_percent}')
            
        if n_less_10 <= n_small_percent:
            return gcad
    return gcad
        

def AdaptiveGroupedAgeCount(cad, bin_widths=[1, 2, 4, 8, 16, 32], max_age=90, low_count=9, lost_percent_limit=0.05, verbose=False):
    """ Adaptively adjust the grouped concept-age distribution bin width 
    
    The idea is to increase the bin widths until the number of patients lost by excluding bins with low counts
    is less than lost_percent_limit% of the total number of patients with the concept
    
    Also, make sure that if patients are lost, at least 5 are lost so that users can't infer any bin having less than 5 patients.
    """
    n_patients = np.sum(cad.counts)
#     n_limit = max(n_patients * lost_percent_limit, 1)
    n_limit = n_patients * lost_percent_limit
    for bw in bin_widths:
        gcad = GroupedAgeCounts(cad, bin_width=bw, max_age=max_age)
        n_lost = np.sum(gcad.counts[gcad.counts <= low_count])
#         n_lost = np.sum((gcad.counts <= low_count) & (gcad.counts > 0)) * low_count
        
        if verbose:
            print(f'bin_width: {bw}\tn_limit: {n_limit}\tn_lost: {n_lost}')
        
#         if n_lost == 0 or (n_lost >= 5 and n_lost <= n_limit):       
        if n_lost <= n_limit:       
            return gcad

    return gcad


def poisson_perturbation(x, n_samples=9):
    """ Applies Poisson perturbation to the input x 
    
    Params
    ------
    x: number, list, or numpy array
    n_samples: For each value in x, take n_samples draws from a Poisson distribution with lambda=x[i] and return the median.
               Note: n_samples should be odd.
    
    Returns
    -------
    values perturbed
    """
    assert n_samples % 2 == 1, 'n_samples must be odd (otherwise creates negative bias)'
    
    try:
        perturbed = [int(np.median(np.random.poisson(xi, n_samples))) for xi in x]

        # If input was numpy array, return an numpy array
        if isinstance(x, np.ndarray):
            perturbed = np.array(perturbed, dtype=np.uint32)
    except TypeError:
        # x is scalar and not iterable
        perturbed = int(np.median(np.random.poisson(x, n_samples)))
        
    return perturbed


class GroupedDeltaCounts():
    """ Groups the DeltaCounter into larger bins """
    def __init__(self, delta_counter, bin_width=2, n=None):
        self.bin_width = bin_width
        self.n = n if n is not None else int(np.ceil(delta_counter.n / bin_width))
        self.bins = self.n * 2 + 1  # Number of bins        
        
        if bin_width == 1 and self.n == delta_counter.n:
            # No change in structure, just change from list to ndarray
            self.counts = np.array(delta_counter.counts)
        else:
            self.counts = np.zeros(self.bins, dtype=np.uint32)
            dcn = delta_counter.n
            cnts = np.array(delta_counter.counts)
            
            # No grouping for 0-day co-occurrence
            self.counts[self.n] = cnts[dcn]
            
            # If the binning stretches "beyond" the original counts array, pad the original counts array
            reach = self.bin_width * self.n
            if reach > dcn:
                pad = np.zeros(reach - dcn, dtype=np.uint32)
                cnts = np.concatenate((pad, cnts, pad))
            dcc = int(np.floor(len(cnts) / 2))
            
            # Fill in the positive deltas            
            upper = dcc + reach + 1
            self.counts[(self.n + 1):(self.bins)] = cnts[(dcc + 1):(upper)].reshape(self.bin_width, self.n, order='F').sum(axis=0)
                        
            # Fill in the negative deltas
            lower = dcc - reach
            self.counts[0:self.n] = cnts[lower:dcc].reshape(self.bin_width, self.n, order='F').sum(axis=0)
            
            # Add the leftover bins
            if reach < dcn:
                self.counts[self.bins - 1] += cnts[upper:]
                self.counts[0] += cnts[:lower]
        
    def bin_labels(self):
        labels = [''] * self.bins
        labels[self.n] = '0'
        
        if self.n > 0:
            labels[slice(self.n - 1, self.n + 2, 2)] = ['-1', '1']
            max_bin_label = f'{2 ** ((self.n-1) * self.bin_width)}+'
            labels[self.bins - 1] = max_bin_label
            labels[0] = '-' + max_bin_label                        
        
        for i in range(2, self.n):
            lower = 2 ** ((i - 1) * self.bin_width)
            upper = 2 ** (i * self.bin_width) - 1
            labels[self.n + i] = f'{lower} — {upper}'
            labels[self.n - i] = f'-{upper} — -{lower}'
            
        return labels
    
    def bin_labels_mixed(self):
        """ Compact labels with mixed units (easier to interpret) """
        labels = ['0d']
        
        for i in range(0, self.n):
            days = 2 ** (i * self.bin_width)
            if days >= 365:
                years = days / 365.25
                label = f'{years:0.1f}y'
            else:
                label = f'{days}d'
            
            labels.append(label)
            labels.insert(0, '-' + label)
            
        return labels
    
    def x(self):
        """ X-ticks ranging from -self.n:self.n"""
        return np.arange(-self.n, self.n+1)
    
    def get(self, ix):
        """ Gets the counts using a relative index, where 0 is the 0-day bin, 1 is the 1-day bin, 
        -1 is the -1-day bin, etc. 
        
        Params
        ------
        indices: relative index or list of indices """        
        try:
            # Assume ix is an iterable and retrieve all requested indices
            if isinstance(ix, list):
                ix = np.array(ix)
            if not isinstance(ix, np.ndarray): 
                raise TypeError('ix expected to be list or numpy.ndarray')
                
            return self.counts[self.n + ix]
        except TypeError:
            # ix is not iterable, assume it's an int
            return self.counts[self.n + ix]
        

def AdaptiveGroupedDeltaCount(delta, settings=[(1, 13), (2, 6), (4, 3), (8, 2), (16, 1)], 
                              low_count=9, n_limit=1, verbose=False):
    """ Adaptively adjust the grouped delta distribution bin width 
    
    The idea is to increase the bin widths until there is at most 1 bin with low_count (excluding bins with 0 count).
    The 0-day bin doesn't contribute to the lost count since it's not included in binning.
    """
    n_patients = np.sum(delta.counts)
    for setting in settings:
        gdc = GroupedDeltaCounts(delta, bin_width=setting[0], n=setting[1])
        low_count_bool_index = (gdc.counts <= low_count) & (gdc.counts > 0)
        low_count_bool_index[gdc.n] = False  # Don't include day-0 bin in the count of lost patients
        
        n_lost = np.sum(low_count_bool_index)
        
        if verbose:
            print(f'settings: {setting}\tn_limit: {n_limit}\tn_lost: {n_lost}')        
        
        if n_lost <= n_limit:
            return gdc
    return gdc        


def AdaptiveGroupedDeltaCount2(delta, settings=[(1, 13), (2, 6), (4, 3), (8, 2), (16, 1)], 
                              low_count=9, lost_percent_limit=0.02, verbose=False):
    """ Adaptively adjust the grouped delta distribution bin width 
    
    The idea is to increase the bin widths until the number of patients lost by excluding bins with low counts
    is less than lost_percent_limit% of the total number of patients with the concept
    """
    n_patients = np.sum(delta.counts)
#     n_limit = max(n_patients * lost_percent_limit, 1)
    n_limit = n_patients * lost_percent_limit
    for setting in settings:
        gdc = GroupedDeltaCounts(delta, bin_width=setting[0], n=setting[1])
        # Don't include day-0 bin in the count of lost patients
        day0_mask = [True] * gdc.n + [False] + [True] * gdc.n
        n_lost = np.sum(gdc.counts[(gdc.counts <= low_count) & day0_mask])
        
        if verbose:
            print(f'settings: {setting}\tn_limit: {n_limit}\tn_lost: {n_lost}')        
        
        if n_lost <= n_limit:
            return gdc
    
    return gdc     