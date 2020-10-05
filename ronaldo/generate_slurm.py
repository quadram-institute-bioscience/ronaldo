#!/usr/bin/env python3
"""
generate_slurm makes submissions scripts, for a dir of results dirs 

generate_slurm makes submissions scripts, for a dir of results dirs 

### CHANGE LOG ### 
2020-09-08 Nabil-Fareed Alikhan <nabil@happykhan.com>
    * Initial build
"""

import time
import os 
import argparse
import logging 
from os import path
import meta
import sys 
import re
epi = "Licence: " + meta.__licence__ +  " by " +meta.__author__ + " <" +meta.__author_email__ + ">"
logging.basicConfig()
log = logging.getLogger()

def main(args):
    log.info(args)
    for dir in os.listdir(args.data_dir):
        illumina_read_path = path.join(args.data_dir, dir, 'ncovIllumina_sequenceAnalysis_trimPrimerSequences/')
        ont_read_path  = path.join(args.data_dir, dir, 'articNcovNanopore_sequenceAnalysisMedaka_articMinIONMedaka')
        valid_read_path = None
        ont = False
        if path.exists(illumina_read_path):
            valid_read_path = illumina_read_path
        elif path.exists(ont_read_path):
            valid_read_path = ont_read_path
            ont = True
        runname_match = re.match('.+\.(20\d+).*', dir)
        if runname_match and valid_read_path:
            runname = runname_match.group(1)
            blank_list = []
            for bam_file in os.listdir(valid_read_path):
                for blank_name in args.blank_prefix:                        
                    if bam_file.lower().startswith(blank_name.lower()) and bam_file.endswith('sorted.bam'):
                        blank_list.append(bam_file)
                        break
            if len(blank_list) > 0: 
                write_pbs(args.output, runname, valid_read_path, blank_list ,args.ctdata, ont, tempdir=args.tempdir)
            else: 
                log.error('no blanks found for ' + runname)

def write_pbs(output_dir, runname, datadir, blanks, ctdata, ont=False, tempdir=None):
    if not path.exists(output_dir):
        os.mkdir(output_dir)
    plat = 'illumina'
    if ont:
        plat = 'nanopore'
    output_script_path  = path.join(output_dir, f'ronaldo.{plat}.{runname}.sh')
    with open(output_script_path, 'w') as out_handle:
        out_handle.write('#!/bin/bash\n')
        out_handle.write('#SBATCH -p qib-long,nbi-short,nbi-medium,nbi-long\n')
        out_handle.write('#SBATCH -t 0-2:00\n')
        out_handle.write('#SBATCH -c 2\n')
        out_handle.write(f'#SBATCH -J ronaldo.{plat}.{runname}\n')
        out_handle.write(f'#SBATCH -J ronaldo.{plat}.{runname}\n')
        out_handle.write('source python-3.7.2\n')
        out_handle.write('source ~/ronaldo/venv/bin/activate\n')
        out_handle.write('cd ~/ronaldo/\n')
        blanks = ' '.join(blanks)
        cmd = 'python ronaldo/ronaldo.py calculate' 
        if tempdir: 
            cmd += f' --tempdir {tempdir} '
        if ont:
            cmd += ' --ont '
        cmd += f'  --ctdata {ctdata} {runname} {datadir} {blanks}'
        out_handle.write(cmd + '\n')

if __name__ == '__main__':
    start_time = time.time()
    log.setLevel(logging.INFO)
    desc = __doc__.split('\n\n')[1].strip()
    parser = argparse.ArgumentParser(description=desc,epilog=epi)
    parser.add_argument ('-v', '--verbose', action='store_true', default=False, help='verbose output')
    parser.add_argument('--version', action='version', version='%(prog)s ' + meta.__version__)
    parser.add_argument('-o','--output',action='store',help='output directory', default='ronaldo_out')
    parser.add_argument('data_dir', action='store', help='Dir with data dirs')
    parser.add_argument('ctdata', action='store', help='ct_data file')
    parser.add_argument('blank_prefix', metavar='N', nargs='+',  help='Prefix to find blanks')
    parser.add_argument('--tempdir',action='store',help='temp directory')

    args = parser.parse_args()
    if args.verbose: 
        log.setLevel(logging.DEBUG)
        log.debug( "Executing @ %s\n"  %time.asctime())    
    main(args)
    if args.verbose: 
        log.debug("Ended @ %s\n"  %time.asctime())
        log.debug('total time in minutes: %d\n' %((time.time() - start_time) / 60.0))
    sys.exit(0)