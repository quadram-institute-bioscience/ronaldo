import csv 
import logging
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter

log = logging.getLogger(__name__)


""" def make_qc_plot(depth_pos, n_density, samplename, window=200):
    depth_df = pd.DataFrame( { 'position' : [pos[1] for pos in depth_pos], 'depth' : [dep[2] for dep in depth_pos] } )
    depth_df['depth_moving_average'] = depth_df.iloc[:,1].rolling(window=window).mean()

    n_df = pd.DataFrame( { 'position' : [pos[0] for pos in n_density], 'n_density' : [dens[1] for dens in n_density] } )

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    ax1.set_xlabel('Position')

    ax1.set_ylabel('Depth', color = 'g')
    ax1.set_ylim(top=10**5, bottom=1)
    ax1.set_yscale('log')
    ax1.plot(depth_df['depth_moving_average'], color = 'g')

    ax2.set_ylabel('N density', color = 'r')  
    ax2.plot(n_df['n_density'], color = 'r')
    ax2.set_ylim(top=1)

    plt.title(samplename)
    plt.savefig(samplename + '.depth.png') """

def fetch_data(filename):
    values = {}
    for x in csv.DictReader(open(filename), dialect=csv.excel):
        values[x['sample_name']] = x 
        values[x['sample_name']]['site'] = filename.split('.')[1]
    return values

def platform_plot(data, output_dir):
    platform_count_1 = dict(Counter([x['ct_platform_1'] for x in data]))
    bp = plt.bar(platform_count_1.keys(), platform_count_1.values())
    platform_count_2 = Counter([x['ct_platform_2'] for x in data])
    bp2 = plt.bar(platform_count_2.keys(), platform_count_2.values())
    plt.xlabel('Diagnostic platform')
    plt.ylabel('Count')    
    plt.savefig(output_dir  + '/ronaldo.platform_plot.png')
    plt.savefig(output_dir  + '/ronaldo.platform_plot.svg')    

def ct_plot(data, output_dir):

    ct = {}
    ct_recovery = {}
    min_value = 25
    max_value = 38
    for x in data:   
        ct_max = round(float(x.get('min_ct_value')))
        if ct_max > 0:
            if ct_max <= min_value:
                ct_max = min_value
            if ct_max >= max_value:
                ct_max = max_value        
            recovery = float(x.get('pc_pos_gte_20'))
            if x.get('sequencing_platform') == 'ILLUMINA':
                recovery = float(x.get('pc_pos_gte_10'))
            if ct.get(ct_max):
                ct[ct_max].append(recovery)
            else:
                ct[ct_max]= [recovery]

    values = [ct[x] for x in sorted(ct)]
    fig, ax = plt.subplots()

    bp = plt.boxplot(values)
    plt.xlabel('CT value')
    plt.ylabel('Genome recovery (%)')
    ticks = [] 
    for x in sorted(ct):
        if x == min_value:
            ticks.append('<' + str(x))
        elif x == max_value:
            ticks.append('>' + str(x))        
        else:    
            ticks.append(str(x))
    ax.set_xticklabels(ticks)
    plt.savefig(output_dir  + '/ronaldo.ct_plot.png')
    plt.savefig(output_dir  + '/ronaldo.ct_plot.svg')
    plt.close()

