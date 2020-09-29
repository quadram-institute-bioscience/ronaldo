#!/usr/bin/env python3
"""
Trash script to locate files for a given sample

Because naming things ten different ways is just what we do.

### CHANGE LOG ### 
2020-09-03 Nabil-Fareed Alikhan <nabil@happykhan.com>
    * Initial build
"""
import time
import argparse
import logging
import meta 
import csv 
import os 
import sys 

epi = "Licence: " + meta.__licence__ +  " by " +meta.__author__ + " <" +meta.__author_email__ + ">"
logging.basicConfig()
log = logging.getLogger()

def main(args): 
    all_samples = {} 
    result_dirs = sorted([x for x in os.listdir(args.datadir) if x.startswith('result.illumina')], reverse=True)
    
    for record in csv.DictReader(open(args.datfile), dialect=csv.excel_tab): 
        cts = [record.get("ct_1_ct_value", 0.0), record.get("ct_2_ct_value", 0.0)]
        new_dict = dict(sequencing_platform = "ILLUMINA", sample_name=record["central_sample_id"], ct_platform_1 = record.get('ct_1_test_platform', 'UNKNOWN'), ct_platform_2 = record.get('ct_2_test_platform', 'UNKNOWN'), max_ct_value=max(cts), min_ct_value=min(cts))
        if new_dict['ct_platform_1'] == '':
            new_dict['ct_platform_1'] = 'UNKNOWN'
        if new_dict['ct_platform_2'] == '':
            new_dict['ct_platform_2'] = 'UNKNOWN'            
        if new_dict['max_ct_value'] == '':
            new_dict['max_ct_value'] = '0'
        if new_dict['min_ct_value'] == '':
            new_dict['min_ct_value'] = '0'
        for dir in result_dirs: 
            file_dir_path =  os.path.join(args.datadir, dir, "ncovIllumina_sequenceAnalysis_readMapping")
            sample_name = record["central_sample_id"].split('-')[1]
            filename = [x for x in os.listdir(file_dir_path) if x.startswith(sample_name)]
            if filename:
                if filename[0] != '':
                    new_dict['filename'] = filename[0]
                    break
        if new_dict.get('filename'):
            all_samples[record["central_sample_id"]] = new_dict
    db_path = 'fixed.dat'
    db_out = csv.DictWriter(open(db_path, 'w'), fieldnames=list(all_samples.values())[0].keys())
    db_out.writeheader()
    db_out.writerows(all_samples.values())



if __name__ == '__main__':
    start_time = time.time()
    log.setLevel(logging.INFO)
    desc = __doc__.split('\n\n')[1].strip()
    parser = argparse.ArgumentParser(description=desc,epilog=epi)
    parser.add_argument ('-v', '--verbose', action='store_true', default=False, help='verbose output')
    parser.add_argument('--version', action='version', version='%(prog)s ' + meta.__version__)
    parser.add_argument('--datfile', action='store', default="temp.dat")
    parser.add_argument('--datadir', action='store', default="/home/ubuntu/transfer/incoming/QIB_Sequencing/Covid-19_Seq/")

    args = parser.parse_args()
    if args.verbose: 
        log.setLevel(logging.DEBUG)
        log.debug( "Executing @ %s\n"  %time.asctime())    
    main(args)
    if args.verbose: 
        log.debug("Ended @ %s\n"  %time.asctime())
        log.debug('total time in minutes: %d\n' %((time.time() - start_time) / 60.0))
    sys.exit(0)