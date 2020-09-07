import pysam
import logging
from collections import Counter

def get_genome_metrics(bam_file):
    recovery_20 = 0 
    recovery_10 = 0 
    total_bases = 0
    coverage = [ ]
    for coord_line in pysam.depth(bam_file ,'-a', '-d', '0').split('\n'):
        coord = coord_line.split('\t')
        if len(coord) > 2:
            total_bases += 1
            coverage.append(float(coord[2]))
            if int(coord[2]) >= 10:
                recovery_10 += 1 
            if int(coord[2]) > 20:
                recovery_20 += 1 
    if total_bases != 29903:
        logging.warn(f'SARSCOV reference genome not expected size, found {total_bases}')
    return round(float(recovery_10) / total_bases * 100, 2) , round(float(recovery_20) / total_bases * 100, 2), round(float(sum(coverage)) / len(coverage),3)


def get_well_mapped_reads(bam_file, read_length=148):
    any_mapped_reads = 0
    bam = pysam.AlignmentFile(bam_file, "rb")
    no_reads = 0
    total_bases = 0 
    for read in bam:
        if not read.is_unmapped:
            any_mapped_reads += 1
            total_bases += len(read.seq)
            match_lengths = [match[1] - match[0] for match in read.cigar]
            if sum(match_lengths) >= read_length:
                no_reads += 1
    return  no_reads, any_mapped_reads, total_bases
