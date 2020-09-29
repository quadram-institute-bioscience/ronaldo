#!/usr/bin/env python3
"""
ronaldo calculate false positive rates for different assays

ronaldo calculate false positive rates for different assays

### CHANGE LOG ### 
2020-09-03 Nabil-Fareed Alikhan <nabil@happykhan.com>
    * Initial build - split from dirty scripts
2020-09-09 Nabil-Fareed Alikhan <nabil@happykhan.com>
    * Fixes for cases where NO reads map. => div0 errors
"""
import collections
import logging
import csv 
import datetime
from os import path, mkdir, listdir
import os 
import argparse
import meta
import sys
import time
from sam_util import get_genome_metrics_nanopore, get_genome_metrics_illumina, get_well_mapped_reads

epi = "Licence: " + meta.__licence__ +  " by " +meta.__author__ + " <" +meta.__author_email__ + ">"
logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s')
log = logging.getLogger(__name__)

def check_blanks(blank_list, ref_length=29000, read_length=148, platform="ILLUMINA"):
    max_coverage = 0
    max_recovery_10 = 0
    max_recovery_20 = 0
    max_reads = 0
    # TODO: Logic to handle Illumina or Nanopore
    for bam_file in blank_list:
        if platform == 'ILLUMINA':
            current_recovery_10, current_recovery_20, current_coverage  = get_genome_metrics_illumina(bam_file)
        else:
            current_recovery_10, current_recovery_20, current_coverage  = get_genome_metrics_nanopore(bam_file)
        current_reads, current_mapped_reads, current_total_bases = get_well_mapped_reads(bam_file, read_length)
        if current_recovery_10 > max_recovery_10:
            max_recovery_10 = current_recovery_10
            max_recovery_20 = current_recovery_20
        if current_coverage > max_coverage:
            max_coverage = current_coverage
        if current_reads > max_reads:
            max_reads = current_reads                        
    return max_coverage, max_recovery_10, max_recovery_20, max_reads

def get_sample_metrics(bam_file, platform='OXFORD_NANOPORE'):
    coverage = 0
    recovery_10 = 0
    recovery_20 = 0
    reads = 0
    recovery_10, recovery_20, coverage = get_genome_metrics(bam_file)
    if platform == 'ILLUMINA':
        reads, any_mapped_reads, total_bases = get_well_mapped_reads(bam_file)
    return coverage, recovery_10, recovery_20, reads

def calculate_metrics(args): 
    # prepopulate ct data from input table 
    log.info(f'Starting RonaLDO on {args.runname}')
    existing_sample_info = {}
    output_sample_info = []
    platform = 'ILLUMINA'
    if args.ont:
        platform = 'OXFORD_NANOPORE'
    if args.ctdata:
        for record in csv.DictReader(open(args.ctdata), dialect=csv.excel): 
            existing_sample_info[record["filename"]] = dict(runname=args.runname, filename=record["filename"], sequencing_platform = platform, sample_name=record["sample_name"], ct_platform_1 = record.get('ct_platform_1', 'UNKNOWN'), ct_platform_2 = record.get('ct_platform_2', 'UNKNOWN'), max_ct_value=record.get("max_ct_value", 0), min_ct_value=record.get('max_ct_value',0))
    else:
        for bam_file in [path.join(args.bamfolder, bam_file) for bam_file in listdir(args.bamfolder) if bam_file not in args.blankbam]: 
            filename = path.basename(bam_file)
            existing_sample_info[filename] = dict(runname=args.runname, sequencing_platform = platform, filename=filename, sample_name=filename, ct_platform_1 = 'UNKNOWN', ct_platform_2 = 'UNKNOWN', max_ct_value=0.0, min_ct_value=0.0)

    # Create path for blanks
    blank_paths = [path.join(args.bamfolder, blank_file) for blank_file in args.blankbam]
    input_good = True
    for blank_path in blank_paths:
        if not path.exists(blank_path):
            log.error(f'Filepath to Blank does not exist: {blank_path}')
            input_good = False
    if input_good:
        #Create output dir 
        if not path.exists(args.db):
            mkdir(args.db)
        # Check blanks
        blank_coverage, blank_recovery_10, blank_recovery_20, blank_reads = check_blanks(blank_paths)
        log.info(f'Max blank genome coverage: {blank_coverage}')
        log.info(f'Max blank genome recovery >10X: {blank_recovery_10}')
        log.info(f'Max blank genome recovery >20X: {blank_recovery_20}')
        log.info(f'Max blank number of mapped reads: {blank_reads}')
        cov_cut_not_ok = False
        if blank_recovery_20 > args.blank_recovery_cutoff and args.ont:
            cov_cut_not_ok = True
        elif blank_recovery_10 > args.blank_recovery_cutoff and not args.ont:
            cov_cut_not_ok = True
        if cov_cut_not_ok or blank_reads > args.blank_read_cutoff:
            log.info('RUN SKIPPED: Specified Blanks from this run have too much SARSCOV2 content')
        else:
            log.info('BLANKS OK!')
            for bam_file in [path.join(args.bamfolder, bam_file) for bam_file in listdir(args.bamfolder) if bam_file not in args.blankbam]: 
                coverage, recovery_10, recovery_20, reads = get_sample_metrics(bam_file, platform=platform)
                bam_filename = path.basename(bam_file)
                new_output_sample_info = existing_sample_info.get(bam_filename)
                if new_output_sample_info:
                    new_output_sample_info.update(dict(blank_coverage=blank_coverage, blank_recovery_10=blank_recovery_10, blank_recovery_20=blank_recovery_20,blank_reads=blank_reads))
                    new_output_sample_info.update(dict(mean_cov=coverage, pc_pos_gte_20=recovery_20, pc_pos_gte_10=recovery_10, no_reads=reads))
                    output_sample_info.append(new_output_sample_info)
                else:
                    log.warning(f'No data for {bam_filename}')
            if output_sample_info:
                db_path = path.join(args.db, f'ronaldo.db.{args.runname}.csv')
                db_out = csv.DictWriter(open(db_path, 'w'), fieldnames=output_sample_info[0].keys())
                db_out.writeheader()
                db_out.writerows(output_sample_info)
            else:
                log.warning('No data found in table')

def assess_run(args):
    # Read dir. 
    all_records = {} 
    for data_file in [x for x in os.listdir(args.db) if x.endswith('.csv')]:
        data_file_path = path.join(args.db, data_file)
        with open(data_file_path) as data_file_handle:
            # pull all tables out and all records. 
            for record in csv.DictReader(data_file_handle, dialect=csv.excel):
                all_records[record['sample_name']] = record
    # Apply cutoffs. 
    for record in all_records.values():
        failed = 0
        record['false_positive'] = False
        # We need to handle where no. blank reads are 0 (esp. for Nanopore) - Suggested by Matt Loose
        corrected_blank_recovery_20 = float(record.get('blank_recovery_20'))
        if corrected_blank_recovery_20 == 0:
            corrected_blank_recovery_20 = 1.00

        if float(record.get('mean_cov')) * args.coverage < float(record.get('blank_coverage')):
            failed += 1        
        if record.get('sequencing_platform') == 'ILLUMINA':
            if float(record.get('pc_pos_gte_10'))  < float(record.get('blank_recovery_10')) * args.recovery:
                failed += 1
            if float(record.get('no_reads')) < float(record.get('blank_reads')) * args.noreads or float(record.get('no_reads')) < args.totalreads:
                failed += 1
            if failed == 3: 
                record['false_positive'] = True
        else:
            if float(record.get('pc_pos_gte_20')) < corrected_blank_recovery_20 * args.recovery:
                failed += 1
            if failed == 2: 
                record['false_positive'] = True

    # Update & output final table with false positives 
    if not path.exists(args.output):
        os.mkdir(args.output)
    out_file_path = path.join(args.output, f'ronaldo.{args.sitename}.summary.csv')
    if all_records:
        with open(out_file_path, 'w') as out_handle:
            out_dict = csv.DictWriter(out_handle, fieldnames=list(all_records.values())[0].keys())
            out_dict.writeheader()
            out_dict.writerows(all_records.values())
    else:
        log.warning('No data found in table')

def is_valid_dir(parser, arg):
    if not path.exists(arg):
        parser.error("The directory %s does not exist!" % arg)
    else:
        return arg


if __name__ == '__main__':
    start_time = time.time()
    log.setLevel(logging.INFO)
    desc = __doc__.split('\n\n')[1].strip()
    parser = argparse.ArgumentParser(description=desc,epilog=epi)
    parser.add_argument ('-v', '--verbose', action='store_true', default=False, help='verbose output')
    parser.add_argument('--version', action='version', version='%(prog)s ' + meta.__version__)
    subparsers = parser.add_subparsers(help='commands')

    # Metrics parser
    metrics_parser = subparsers.add_parser('calculate', help='Calculate genome cov metrics')
    metrics_parser.add_argument('-d','--db', action='store',help='DB directory', default='ronaldo_db')
    metrics_parser.add_argument('--ctdata', action='store',  help='Path to table with assay information')
    metrics_parser.add_argument('--blank_read_cutoff', action='store',  default=500,  help='Run skipped if blanks have number of mapped reads')
    metrics_parser.add_argument('--blank_recovery_cutoff', action='store', default=4.0, help='Run skipped if blanks have higher perc. genome recovery')
    metrics_parser.add_argument('--ont', action='store_true', default=False, help='Data is OXFORD NANOPORE')
    metrics_parser.add_argument('-l','--readlen',action='store',help='Minimum length for a mapped read (use with ILLUMINA only)', default=148)
    metrics_parser.add_argument('runname', action='store', help='Informative label for this run')    
    metrics_parser.add_argument('bamfolder', action='store', help='Folder of SARSCOV2 BAM files', type=lambda x: is_valid_dir(metrics_parser, x))
    metrics_parser.add_argument('blankbam', metavar='N', nargs='+', help='Negative control BAM file')
    metrics_parser.set_defaults(func=calculate_metrics)

    # Filter parser
    filter_parser = subparsers.add_parser('filter', help='Filter metric results, determine false positives')
    filter_parser.add_argument('-d','--db', action='store',help='DB directory', default='ronaldo_db')
    filter_parser.add_argument('-o','--output',action='store',help='output directory', default='ronaldo_out')
    filter_parser.add_argument('-c','--coverage',action='store',help='Minimum fold genome coverage for mapped reads', default=2)
    filter_parser.add_argument('-r','--recovery',action='store',help='Minimum fold genome recovery for mapped reads', default=2)    
    filter_parser.add_argument('-n','--noreads',action='store',help='Minimum fold number of mapped reads', default=5)
    filter_parser.add_argument('-t','--totalreads',action='store',help='Minimum total number of mapped reads', default=30)
    filter_parser.add_argument('sitename', action='store', help='Informative label for your site')    
    filter_parser.set_defaults(func=assess_run)

    args = parser.parse_args()
    if args.verbose: 
        log.setLevel(logging.DEBUG)
        log.debug( "Executing @ %s\n"  %time.asctime())    
    if hasattr(args, 'func'):
        args.func(args)
    else: 
        parser.print_help()
    if args.verbose: 
        log.debug("Ended @ %s\n"  %time.asctime())
        log.debug('total time in minutes: %d\n' %((time.time() - start_time) / 60.0))
    sys.exit(0)