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
    blacklist = ['result.illumina.20200930.depricated', 'result.ont.cog35.200929.single_end']
    for record in csv.DictReader(open(args.datfile), dialect=csv.excel_tab): 
        cts = [record.get("ct_1_ct_value", 0.0), record.get("ct_2_ct_value", 0.0)]
        new_dict = dict(sequencing_platform=record['instrument_make'], sample_name=record["central_sample_id"], ct_platform_1 = record.get('ct_1_test_platform', 'UNKNOWN'), ct_platform_2 = record.get('ct_2_test_platform', 'UNKNOWN'), max_ct_value=max(cts), min_ct_value=min(cts))
        if new_dict['ct_platform_1'] == '':
            new_dict['ct_platform_1'] = 'UNKNOWN'
        if new_dict['ct_platform_2'] == '':
            new_dict['ct_platform_2'] = 'UNKNOWN'            
        if new_dict['max_ct_value'] == '':
            new_dict['max_ct_value'] = '0'
        if new_dict['min_ct_value'] == '':
            new_dict['min_ct_value'] = '0'
        all_samples[record["central_sample_id"]] = new_dict
    result_dirs = sorted([x for x in os.listdir(args.datadir) if x.startswith('result.') and x not in blacklist ], reverse=True)
    for dir in result_dirs: 
        data_dir  =  os.path.join(args.datadir, dir)
        qc_file = [x for x in os.listdir(data_dir) if x.endswith('qc.csv')]
        if qc_file:
            samples = [x for x in csv.DictReader(open(os.path.join(args.datadir, dir, qc_file[0]))) ] 
            if dir.startswith('result.illumina'):
                seq_plat = 'ILLUMINA'
                file_dir_path =  os.path.join(args.datadir, dir, "ncovIllumina_sequenceAnalysis_trimPrimerSequences")
            else:
                seq_plat = 'OXFORD_NANOPORE'      
                file_dir_path =  os.path.join(args.datadir, dir, "articNcovNanopore_sequenceAnalysisMedaka_articMinIONMedaka")            
            for sample in samples:
                if seq_plat == 'OXFORD_NANOPORE':
                    sample_name = sample['sample_name']
                else:
                    sample_name = sample['sample_name'].split('_')[0]
                    if not sample_name.startswith('NORW'):
                        sample_name = 'NORW-' + sample_name                    
                if all_samples.get(sample_name):
                    bam_file_path = os.path.join(file_dir_path, sample['bam'])
                    if os.path.exists(bam_file_path):
                        all_samples[sample_name]['filename'] = sample['bam']
                    else: 
                        log.info('Path not found ' + bam_file_path)

        else:
            log.info(f'No QC file in {dir}')
    db_path = 'fixed.dat'
    good_files = [x for x in all_samples.values() if x.get('filename')]
    db_out = csv.DictWriter(open(db_path, 'w'), fieldnames=good_files[0].keys())
    db_out.writeheader()
    db_out.writerows(good_files)


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