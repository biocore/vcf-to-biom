#!/usr/bin/env python

from numpy import inf
from sys import argv
from biom.util import biom_open
from biom.table import table_factory

def biom_data_from_vcfs(vcfs, min_position=0, max_position=inf):
    oids = {}
    ordered_oids = []
    sids = {}
    ordered_sids = []
    data = {}
    master_oids = set([])
    for vcf in vcfs:
        working_oids = set([])
        vcf = biom_open(vcf)
        for line in vcf:
            fields = line.strip().split('\t')
            if fields[0] == '#CHROM':
                pass
            elif fields[0].startswith("#"):
                pass
            else:
                chrom = fields[0]
                pos = int(fields[1])
                oid = '%s.%d' % (chrom, pos)
                working_oids.add(oid)
        if len(master_oids) == 0:
            master_oids = working_oids
        else:
            master_oids = set.intersection(master_oids, working_oids)
#            master_oids = master_oids | working_oids
        vcf.close()
    for vcf in vcfs:
        vcf = biom_open(vcf)
        for line in vcf:
            fields = line.strip().split('\t')
            if fields[0] == '#CHROM':
                # this will differ for human data (when multiple genomes per vcf):
                sid = fields[9]
                try:
                    sid_index = sids[sid]
                except KeyError:
                    ordered_sids.append(sid)
                    sid_index = len(ordered_sids) - 1
                    sids[sid] = sid_index
            elif fields[0].startswith("#"):
                pass
            else:
                chrom = fields[0]
                pos = int(fields[1])
                oid = '%s.%d' % (chrom, pos)
                if fields[4] != '.' and \
                   min_position <= pos <= max_position and \
                   oid in master_oids:
                    try:
                       oid_index = oids[oid]
                    except KeyError:
                        ordered_oids.append(oid)
                        oid_index = len(ordered_oids) - 1
                        oids[oid] = oid_index
                    # this will differ for non-haploid data:
                    data[(oid_index, sid_index)] = 1

    return data, ordered_oids, ordered_sids
    
    
if __name__ == "__main__":
    usage = "simple_vcf_to_biom.py output-biom-fp input-vcf-fp1 [input-vcf-fp2 [...]]"
    if len(argv) < 3 or argv[1] == '-h':
        print usage
        exit(0)
    
    output_biom_fp = argv[1]
    vcfs = argv[2:]
    
    data, oids, sids = biom_data_from_vcfs(vcfs)
    print oids, data
    table = table_factory(data, sample_ids=sids, observation_ids=oids)
    
    output_biom_f = open(output_biom_fp,'w')
    table.getBiomFormatJsonString("simple_vcf_to_biom", direct_io=output_biom_f)
    output_biom_f.close()
    
    
    
# data = {(0, 0): 1, (0, 1): 1, (1, 0): 1, (1, 1):0}
# oids = ['10.89673612', '10.89673554']
# #sids = [sample1, sample2]
# {'observation':{}, 'observation2', 0}
# 

############## End Greg's Code################
# 
# {'observation':{'sample1':1, 'sampple2':0}, 'observation2':'sample1':0} 
# 
# 
# def create_biom_table(f, data=None):
#     '''the input is an open vcf file. The output is all of the objects required for
#     table factory.'''
#     if data is None:
#         observation_ids =  = {}
#         sample_ids = {}
#     for line in f:
#         if line.startswith('##'):
#             pass
#         elif line.startswith('INFO'):
#             pass
#         elif line.startswith('#CHROM'): 
#             sample_ids = process_header_entry(line)
#         else:
#             if sample_ids == None:
#                 raise ValueError, "Didn't find '#CHROM' line before data lines.\
#                  Can't continue."
#             else: 
#                 chr_n, pos_n, rs_n, ref, alt, indiv_ids = \
#                 process_data_entry_line(line, sample_ids)
#                 # If length of any snp is greater than one it will be ignored.
#                 if len(ref) > 1 or max(len(x) for x in alt.split(',')) > 1 or len(alt) > 1:
#                     pass
#                 else: 
#                     observation_id = "%s:%s" %(chr_n, pos_n)
#                     for indiv, variation in indiv_ids:
#                         if len(variation) == 2:
#                             if variation == [0, 0]:
#                                 data_row = 0
#                             elif variation == [0, 1]:
#                                 data_row = 1
#                             elif variation == [1, 0]: 
#                                 data_row = 1                        
#                             elif variation == [1, 1]:
#                                 data_row = 2
#                         else:
#                             data_row.append(variation[0])
#                         data[observation_id][indiv] = 
#     if len(data) == 0:
#         raise ValueError, """No valid SNP data was present in the file. Indels will be 
# ignored"""
#     else:
#         data = array(data)
#     return data, sample_ids, observation_ids, sample_md, observation_md





#time python /home/johnchase/vcf_bacteria_data/alpha_code/simple_vcf_to_biom.py /home/johnchase/vcf_bacteria_data/vcfs/test_out.biom /home/johnchase/vcf_bacteria_data/vcfs/Bpseudo-INT2-217_CAGATC_L008_001-novo-gatk.vcf /home/johnchase/vcf_bacteria_data/vcfs/Bpseudo-INT2-87_TTAGGC_L008_001-novo-gatk.vcf /home/johnchase/vcf_bacteria_data/vcfs/Bpseudo-INT2-223_ACTTGA_L008_001-novo-gatk.vcf /home/johnchase/vcf_bacteria_data/vcfs/Bpseudo-INT2-91_TGACCA_L008_001-novo-gatk.vcf /home/johnchase/vcf_bacteria_data/vcfs/Bpseudo-INT2-24_ATCACG_L008_001-novo-gatk.vcf /home/johnchase/vcf_bacteria_data/vcfs/Burkholderia-pseudomallei-INT37Bp018_TTCCATTG_L004_001-novo-gatk.vcf /home/johnchase/vcf_bacteria_data/vcfs/Bpseudo-INT2-38_CGATGT_L008_001-novo-gatk.vcf