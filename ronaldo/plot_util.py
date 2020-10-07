import csv 
import logging
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter 
import os 

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
    #bp2 = plt.bar(platform_count_2.keys(), platform_count_2.values())
    plt.xlabel('Diagnostic platform')
    plt.xticks(rotation=45)
    plt.legend()
    plt.ylabel('Count')    
    plt.savefig(output_dir  + '/ronaldo.platform_plot.png', bbox_inches='tight')
    plt.savefig(output_dir  + '/ronaldo.platform_plot.svg', bbox_inches='tight')
    plt.close()

def platform_fail_plot(data, output_dir, plat_cut=50):
    # Absolute failure count
    fail_platform_count_1 = dict(Counter([x['ct_platform_1'] for x in data if x['false_positive'] == 'True' and x['ct_platform_1'] != 'UNKNOWN']))
    bp = plt.bar(fail_platform_count_1.keys(), fail_platform_count_1.values())
    fail_platform_count_2 = dict(Counter([x['ct_platform_2'] for x in data if x['false_positive'] == 'True']))
    #bp2 = plt.bar(fail_platform_count_2.keys(), fail_platform_count_2.values())
    plt.xlabel('Diagnostic platform')
    plt.xticks(rotation=45)
    plt.legend()
    plt.ylabel('Count')    
    plt.savefig(output_dir  + '/ronaldo.platform_fail_plot.png', bbox_inches='tight')
    plt.savefig(output_dir  + '/ronaldo.platform_fail_plot.svg', bbox_inches='tight')    
    plt.close()
    # pct of total
    fail_dict = []
    platform_count_1 = dict(Counter([x['ct_platform_1'] for x in data]))
    platform_count_2 = Counter([x['ct_platform_2'] for x in data])
    fail_platform_count_1 = dict(Counter([x['ct_platform_1'] for x in data if x['false_positive'] == 'True' and x['ct_platform_1'] != 'UNKNOWN']))
    fail_platform_pct_1 = {}
    for plat, value in fail_platform_count_1.items():
        if platform_count_1[plat] >= plat_cut:
            fail_platform_pct_1[plat] = value / platform_count_1[plat] *100.00
        fail_dict.append(dict(platform=plat,fail_count=value, total_count=platform_count_1[plat] ))
    out_file = os.path.join(output_dir, 'ronaldo.fail_table_plt1.csv')
    tab = csv.DictWriter(open(out_file, 'w'), fieldnames=['platform', 'fail_count', "total_count"])        
    tab.writeheader()
    tab.writerows(fail_dict)    
    bp = plt.bar(fail_platform_pct_1.keys(), fail_platform_pct_1.values())
    fail_dict = []
    fail_platform_count_2 = dict(Counter([x['ct_platform_2'] for x in data if x['false_positive'] == 'True']))
    for plat, value in fail_platform_count_2.items():
        fail_platform_count_2[plat] = value / platform_count_2[plat] *100.00    
        fail_dict.append(dict(platform=plat,fail_count=value, total_count=platform_count_1[plat] ))
    # bp2 = plt.bar(fail_platform_count_2.keys(), fail_platform_count_2.values())
    plt.xlabel('Diagnostic platform')
    plt.xticks(rotation=45)
    plt.legend()
    plt.ylabel('Failure rate of total samples (%)')    
    plt.savefig(output_dir  + '/ronaldo.platform_pct_fail_plot.png', bbox_inches='tight')
    plt.savefig(output_dir  + '/ronaldo.platform_pct_fail_plot.svg', bbox_inches='tight')    
    plt.close()
    out_file = os.path.join(output_dir, 'ronaldo.fail_table_plt2.csv')
    tab = csv.DictWriter(open(out_file, 'w'), fieldnames=['platform', 'fail_count', "total_count"])
    tab.writeheader()
    tab.writerows(fail_dict)


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

