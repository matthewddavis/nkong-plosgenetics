#!/usr/bin/python
import sys
import pybedtools as pbt
import Bio.SeqIO as SeqIO
from Bio.Alphabet.IUPAC import unambiguous_dna
import Bio
import re

# get genome file
genome_fasta_fn = '..//reference/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa'
genome = Bio.SeqIO.to_dict(Bio.SeqIO.parse(genome_fasta_fn, format='fasta', alphabet=unambiguous_dna))

experiments = []
experiments.append({ 'bed_fn' : '../data/reads/chip_exo/pre_b/mef2_goat_preB_chip_september2012/macs/Sample_NKRT_3_index1_R1.trimmed.macs_peaks.bed',
                     'chip_sample_name' : 'mef2_preB',
                 })

experiments.append({ 'bed_fn' : '../data/reads/chip_exo/hpc/mef2_goat_hpc_chip_july2013/macs/Sample_NKRT5A_index2_R1.trimmed.macs_peaks.bed',
                     'chip_sample_name' : 'mef2_hpc',
                 })

experiments.append({ 'bed_fn' : '../data/reads/chip_exo/pre_b/ebf1_goat_preB_chip_february2013/macs/Sample_NKRT4_index3_R1.trimmed.macs_peaks.bed',
                     'chip_sample_name' : 'ebf1_preb',
                 })

genes_gtf_fn = '../genes.gtf'

genes_bt = pbt.BedTool(genes_gtf_fn)

for experiment in experiments:
    chip_sample_name = experiment['chip_sample_name']
    peaks_bed_fn = experiment['bed_fn']

    # get closest peaks
    closest_peaks = genes_bt.closest(peaks_bed_fn, D=True, s=False)

    experiment['closest_peaks'] = closest_peaks
    
    closest_peaks.saveas(experiment['chip_sample_name'] + '.closest_genes.bed')

    up_target_gene_names = ['Ebf1', 'Myb', 'Il7r', 'Ets1', 'Foxo1']
    down_target_gene_names = ['Ciita', 'Fcrl1', 'Il21r', 'Csf2r', 'Flt3l', 'Fcgr3', 'Stat1']

    target_gene_names = []
    target_gene_names.extend(up_target_gene_names)
    target_gene_names.extend(down_target_gene_names)


    coords_list = []

    for gene in closest_peaks:
        if gene.name in target_gene_names:
            coords = dict(zip(['peak_id', 'gene_id', 'chr', 'start', 'stop'], [ gene.fields[12], gene.name, gene.fields[0], int(gene.fields[10]), int(gene.fields[11]) ]))
            if coords not in coords_list:
                coords_list.append(coords)

    for coords in coords_list:
        seq = genome[ coords['chr']][ coords['start'] : coords['stop'] ].seq

            # mef2c pattern
            #pattern = '[AT]{6}TAG'
            # Ebf1 pattern
            #pattern = 'CCC[ACGTN]{2}GGG'
            # test pattern
            #pattern = 'CC[ACTGN]{2}G'
        patterns = [ '[AT]{6}TAG', 'CCC[ACGTN]{2}GGG' ]
        for pattern in patterns:
            hits = re.finditer(pattern, seq.tostring().upper())

            for hit in hits:
                hit_start = hit.start()
                hit_end = hit.end()

                hit_seq = seq[ hit_start : hit_end ].tostring().upper()

                hit_start = hit.start() + coords['start']
                hit_end = hit.end() + coords['start']

                print chip_sample_name, pattern, coords['gene_id'], coords['chr'] + ':' + str(hit_start), '+', hit_seq, "%(peak_id)s %(chr)s:%(start)s-%(stop)s" % (coords)

        rev_seq = seq.reverse_complement()

        for pattern in patterns:
            hits = re.finditer(pattern, rev_seq.tostring().upper())

            for hit in hits:
                hit_start = hit.start()
                hit_end = hit.end()

                hit_seq = rev_seq[ hit_start : hit_end ].tostring().upper()

                hit_end = (len(rev_seq) - hit.start()) +coords['start']
                hit_start = (len(rev_seq) - hit.end()) + coords['start']

                print chip_sample_name, pattern, coords['gene_id'], coords['chr'] + ':' + str(hit_start), '-', hit_seq, "%(peak_id)s %(chr)s:%(start)s-%(stop)s" % (coords)




    
