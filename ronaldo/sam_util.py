import pysam
import logging

def get_genome_recovery(bam_file):
    recovery_20 = 0 
    recovery_10 = 0 
    total_bases = 0
    for coord_line in pysam.depth(bam_file ,'-a', '-d', '0').split('\n'):
        coord = coord_line.split('\t')
        if len(coord) > 2:
            total_bases += 1
            if int(coord[2]) >= 10:
                recovery_10 += 1 
            if int(coord[2]) > 20:
                recovery_20 += 1 
    if total_bases != 29903:
        logging.warn(f'SARSCOV reference genome not expected size, found {total_bases}')
    return round(float(recovery_10) / total_bases * 100, 2) , round(float(recovery_20) / total_bases * 100, 2)


def get_genome_coverage(bam_file):
    no_reads = 0 
    coverage = 0.0
    return coverage, no_reads 
