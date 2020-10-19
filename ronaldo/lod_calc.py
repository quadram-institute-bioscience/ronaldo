import os
import csv 
from sam_util import get_genome_metrics      
import matplotlib.pyplot as plt
from  random import randint
import logging 

logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

LOD_LOC = {'Replicate 3': '/home/ubuntu/transfer/incoming/QIB_Sequencing/Covid-19_Seq/result.illumina.20200819',
    'Replicate 1': "/home/ubuntu/transfer/incoming/QIB_Sequencing/Covid-19_Seq/result.illumina.20200729.LOD1", 
    'Replicate 2': "/home/ubuntu/transfer/incoming/QIB_Sequencing/Covid-19_Seq/result.illumina.20200729.LOD2",
}

all_dat = {}

for name, path  in LOD_LOC.items(): 
    #  Find QC sheet 
    qc_sheet = [os.path.join(path, x) for x in os.listdir(path) if x.endswith('qc.csv')]
    if qc_sheet:
        lod_samples = [x for x in csv.DictReader(open(qc_sheet[0]), dialect=csv.excel) if x['sample_name'].upper().startswith('LOD')]
        # Get file names - make dict of samples 
        for sample in lod_samples:
            if name == 'Replicate 1':
                conc = int(sample['sample_name'].split('_')[0].replace('LoD', ''))
            if name == 'Replicate 2':
                conc = int(sample['sample_name'].split('_')[0].split('-')[2])
            if name == 'Replicate 3':                 
                conc = int(sample['sample_name'].split('_')[0].split('-')[1])
            bam_path = os.path.join(path, "ncovIllumina_sequenceAnalysis_trimPrimerSequences", sample['bam'])
            # Calculate genome recovery & coverage 
            recovery_10, recovery_20, coverage, reads = get_genome_metrics(bam_path)  # (randint(1,100), 1,1,1) #
            if not all_dat.get(name):
                all_dat[name] = {} 
            if not all_dat.get(name).get(conc):
                all_dat[name][conc] = {}
            all_dat[name][conc] = dict(recovery_10=recovery_10, coverage=coverage, reads=reads)
    else:
        log.error('NO QC SHEET ' + qc_sheet)
# plot 
print(all_dat)
plt.style.use('ggplot')
for name, lod in all_dat.items():
    rvalues = [y['recovery_10'] for x,y  in lod.items()]
    lvalues = [x for x,y  in lod.items()]
    plt.scatter(lvalues, rvalues, label=name)

plt.legend()
plt.xlabel('Viral particles per sample (log)')
plt.ylabel('Genome recovery (%)') 
plt.savefig('ronaldo.lod.png', bbox_inches='tight')
plt.savefig('ronaldo.lod.svg', bbox_inches='tight')
plt.close()

import math 

plt.style.use('ggplot')
for name, lod in all_dat.items():
    rvalues = [y['coverage'] for x,y  in lod.items()]
    lvalues = [x for x,y  in lod.items()]
    plt.scatter(lvalues, rvalues, label=name)

plt.legend()
plt.xlabel('Viral particles per sample (log)')
plt.ylabel('Genome coverage (X)') 
plt.savefig('ronaldo.cov_lod.png', bbox_inches='tight')
plt.savefig('ronaldo.cov_lod.svg', bbox_inches='tight')
plt.close()


plt.style.use('ggplot')

for name, lod in all_dat.items():

    rvalues = [y['recovery_10'] for x,y  in lod.items()  ]
    lvalues = [x+1 for x,y  in lod.items() ]
    plt.scatter(lvalues, rvalues, label=name)

plt.legend()
plt.xlabel('Viral particles per sample (log)')
plt.ylabel('Genome recovery (%)') 
plt.xscale('log')
plt.xlim(0, 110)
plt.savefig('ronaldo.lod.log.png', bbox_inches='tight')
plt.savefig('ronaldo.lod.log.svg', bbox_inches='tight')
plt.close()

for name, lod in all_dat.items():
    rvalues = [y['coverage'] for x,y  in lod.items()  ]
    lvalues = [x+1 for x,y  in lod.items() ]
    plt.scatter(lvalues, rvalues, label=name)

plt.legend()
plt.xscale('log')
plt.xlabel('Viral particles per sample (log)')
plt.ylabel('Genome coverage (X)') 
plt.savefig('ronaldo.cov_lod.log.png', bbox_inches='tight')
plt.savefig('ronaldo.cov_lod.log.svg', bbox_inches='tight')
plt.close()