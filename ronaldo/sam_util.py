import pysam
import logging
from collections import Counter

def get_genome_metrics_nanopore(bam_file):
    recovery_20 = 0 
    recovery_10 = 0 
    total_bases = 0
    coverage = [ ]
    try:
        for coord_line in pysam.depth(bam_file ,'-a', '-d', '0').split('\n'):
            coord = coord_line.split('\t')
            if len(coord) > 2:
                total_bases += 1
                coverage.append(float(coord[2]))
                if int(coord[2]) >= 10:
                    recovery_10 += 1 
                if int(coord[2]) > 20:
                    recovery_20 += 1 
    except (pysam.utils.SamtoolsError, OSError) as ex:
        logging.error(f'Error opening BAM file {bam_file}')
        logging.error(ex)
        return 0,0,0
    if total_bases != 29903 and total_bases > 0 :
        logging.warn(f'SARSCOV reference genome not expected size, found {total_bases}')
    if total_bases == 0 or not coverage:
        logging.warn(f'No reads mapped at all in  {bam_file}')
        return 0, 0, 0
    return round(float(recovery_10) / total_bases * 100, 2) , round(float(recovery_20) / total_bases * 100, 2), round(float(sum(coverage)) / len(coverage),3)


def get_genome_metrics_illumina(bam_file):
    recovery_20 = 0 
    recovery_10 = 0 
    total_bases = 0
    coverage = [ ]
    try:
        for coord_line in pysam.depth(bam_file ,'-a', '-d', '0').split('\n'):
            coord = coord_line.split('\t')
            if len(coord) > 2:
                total_bases += 1
                coverage.append(float(coord[2]))
                if int(coord[2]) >= 10:
                    recovery_10 += 1 
                if int(coord[2]) > 20:
                    recovery_20 += 1 
    except (pysam.utils.SamtoolsError, OSError) as ex:
        logging.error(f'Error opening BAM file {bam_file}')
        logging.error(ex)
        return 0,0,0
    if total_bases != 29903 and total_bases > 0 :
        logging.warn(f'SARSCOV reference genome not expected size, found {total_bases}')
    if total_bases == 0 or not coverage:
        logging.warn(f'No reads mapped at all in  {bam_file}')
        return 0, 0, 0
    return round(float(recovery_10) / total_bases * 100, 2) , round(float(recovery_20) / total_bases * 100, 2), round(float(sum(coverage)) / len(coverage),3)


def get_well_mapped_reads(bam_file, read_length=148):
    any_mapped_reads = 0
    no_reads = 0
    total_bases = 0 
    try:    
        bam = pysam.AlignmentFile(bam_file, "rb")
        for read in bam:
            if not read.is_unmapped:
                any_mapped_reads += 1
                total_bases += len(read.seq)
                match_lengths = [match[1] - match[0] for match in read.cigar]
                if sum(match_lengths) >= read_length:
                    no_reads += 1
    except (pysam.utils.SamtoolsError, OSError) as ex:
        logging.error(f'Error opening BAM file {bam_file}')
        logging.error(ex)
        return 0,0,0
    return  no_reads, any_mapped_reads, total_bases

def read_is_fully_mapped(read):

    return True

def get_full_mapped_coverage(bam_file, read_length):
    bam = pysam.AlignmentFile(bam_file, "rb")
    bam.count_coverage(read_callback=read_is_fully_mapped)
