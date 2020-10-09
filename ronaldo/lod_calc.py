import os
import csv 
from sam_util import get_genome_metrics      

import logging 

logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s')
log = logging.getLogger(__name__)

LOD_LOC = {'lod3': '/home/ubuntu/transfer/incoming/QIB_Sequencing/Covid-19_Seq/result.illumina.20200819',
    'lod1': "/home/ubuntu/transfer/incoming/QIB_Sequencing/Covid-19_Seq/result.illumina.20200729.LOD1", 
    'lod2': "/home/ubuntu/transfer/incoming/QIB_Sequencing/Covid-19_Seq/result.illumina.20200729.LOD2",
}

all_dat = {}

for name, path  in LOD_LOC.items(): 
    #  Find QC sheet 
    qc_sheet = [os.path.join(path, x) for x in os.listdir(path) if x.endswith('qc.csv')]
    if qc_sheet:
        lod_samples = [x for x in csv.DictReader(open(qc_sheet[0]), dialect=csv.excel) if x['sample_name'].startswith('LOD')]
        # Get file names - make dict of samples 
        for sample in lod_samples:
            if name == 'lod1':
                conc = sample['sample_name'].split('_')[1].replace('LoD', '')
            if name == 'lod2':
                 conc = sample['sample_name'].split('-')[1]
            bam_path = os.path.join(path, "ncovIllumina_sequenceAnalysis_trimPrimerSequences", sample['bam'])
            # Calculate genome recovery & coverage 
            recovery_10, recovery_20, coverage, reads = get_genome_metrics(bam_path)
            if not all_dat.get(replicate):
                all_dat[replicate] = {} 
            if not all_dat.get(replicate).get(conc):
                all_dat[replicate][conc] = {}
            all_dat[name][conc] = dict(recovery_10=recovery_10, coverage=coverage, reads=reads)
            
# plot 
print(all_dat)

