import pysam
import logging
from collections import Counter
import tempfile
import os 

log = logging.getLogger(__name__)


def get_genome_metrics(bam_file, ref_length = 29903, platform = 'ILLUMINA', read_length=148, verbose=False, temp=None):
    if verbose:
        log.setLevel(logging.DEBUG)
    recovery_20 = 0 
    recovery_10 = 0 
    total_bases = 0
    reads = 0 
    coverage = [ ]
    unmapped_reads = 0 
    fd = None 
    try:
        # If illumina create temp bam of full mapped reads 
        if platform == 'ILLUMINA':
            log.debug('Creating temp file and filtering for ' + os.path.basename(bam_file) )
            if temp:
                fd = tempfile.NamedTemporaryFile(prefix='ronaldo', dir=temp)
            else:
                fd = tempfile.NamedTemporaryFile(prefix='ronaldo')
            bam = pysam.AlignmentFile(bam_file, "rb")
            if not bam.has_index():
                pysam.index(bam_file)
                bam = pysam.AlignmentFile(bam_file, "rb")
            # Pre filter full mapped reads .
            temp_bam = pysam.AlignmentFile(fd, 'wb', template=bam)
            for read in bam.fetch():
                if not read.is_unmapped:
                    match_lengths = [match[1] - match[0] for match in read.cigar]
                    if sum(match_lengths) >= read_length: 
                        temp_bam.write(read)
            bam.close()
            temp_bam.close()
            clean_bam_file = fd.name
        else:
            clean_bam_file = bam_file
        log.debug('Fetching read stats for ' + os.path.basename(bam_file) )
        log.debug('Path: ' + clean_bam_file)
        clean_bam = pysam.AlignmentFile(clean_bam_file, "rb")
        if not clean_bam.has_index(): 
            pysam.index(clean_bam_file)
        for stat_line in pysam.idxstats(clean_bam_file).split('\n'):
            if stat_line.startswith('*'):
                unmapped_reads = int(stat_line.split('\t')[2]) + int(stat_line.split('\t')[3])
            elif stat_line.startswith('MN908947.3'):
                reads = int(stat_line.split('\t')[2])
            else:
                if len(stat_line) > 1:
                    ref_name = stat_line.split('\t')[0]
                    logging.warn(f'Other reference file detected, {ref_name}.')       
        log.debug('Calculating depth for ' + os.path.basename(bam_file) )       
        try:              
            for coord_line in pysam.depth(clean_bam_file ,'-a', '-d', '0').split('\n'):
                coord = coord_line.split('\t')
                if len(coord) > 2:
                    total_bases += 1
                    coverage.append(float(coord[2]))
                    if int(coord[2]) >= 10:
                        recovery_10 += 1 
                    if int(coord[2]) > 20:
                        recovery_20 += 1 
        except (pysam.utils.SamtoolsError, OSError) as ex:
            log.error(f'Error opening BAM file {clean_bam_file}')
            log.error(ex)
            return 0,0,0,0                                    
    except (pysam.utils.SamtoolsError, OSError) as ex:
        log.error(f'Error opening BAM file {bam_file}')
        log.error(ex)
        return 0,0,0,0     
    if platform == 'ILLUMINA' and fd:                    
        fd.close()          
    if total_bases != ref_length and total_bases > 0 :
        log.warn(f'SARSCOV reference genome not expected size, found {total_bases}')
    if total_bases == 0 or not coverage:
        log.warn(f'No reads mapped at all in  {bam_file}')
        return 0, 0, 0, 0
    return round(float(recovery_10) / total_bases * 100, 2) , round(float(recovery_20) / total_bases * 100, 2), round(float(sum(coverage)) / len(coverage),3), reads