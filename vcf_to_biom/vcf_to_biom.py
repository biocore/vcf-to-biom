#!/usr/bin/env python
# File created on 11 Jul 2013
from __future__ import division

__author__ = "John Chase"
__copyright__ = "Copyright 2013"
__credits__ = ["John Chase"]
__license__ = "GPL"
__version__ = "0.0.1-dev"
__maintainer__ = "John Chase"
__email__ = "chasejohnh@gmail.com"
__status__ = "Development"

from biom.table import table_factory, SparseOTUTable
from biom.parse import MetadataMap, generatedby
from numpy import array
from os.path import join
import gzip
from numpy import inf


def process_data_entry_line(line, ids):
    fields = line.strip().split('\t')
    chr_n = fields[0] 
    try:
        pos_n = fields[1]
    except:
        print fields
    rs_n = fields[2]
    ref = fields[3] 
    alt = fields[4]
    indiv_ids = zip(ids, indiv_snp_variation(line))
    return (chr_n, pos_n, rs_n, ref, alt, indiv_ids)

#process the information line. This line contains all of the information about the 
#individuals. The function returns a list of the individuals 
def process_header_entry(line):
    fields = line.split() 
    return fields [9:]

#Creates a list of allelic variation for each SNP. The list itself contains one entry of 
#two numbers for each individual. 
def indiv_snp_variation(line): 
    genotypes = line.split('\t')[9:]
    formatted_genotypes = [] 
    for i in genotypes:
        try: 
            formatted_genotypes.append(map(int,i.split(':')[0].split('|')))
        except:
            try:
                formatted_genotypes.append(map(int,i.split(':')[0].split('/')))
            except:
                formatted_genotypes.append([0])
    return formatted_genotypes
 
# def create_biom_table(f):
#     '''the input is an open vcf file. The output is all of the objects required for
#     table factory.'''
#     sample_ids = None
#     observation_ids = []
#     observation_md = []
#     data = []
#     sample_md = None
#     previous_id = None
#     c =0
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
#                     if observation_id == previous_id:
#                         pass
#                     elif alt == '.':
#                         pass
#                     else:   
#                         observation_ids.append(observation_id)
#                         meta_dic = {"alleles":(ref, alt),"rs":rs_n}
#                         observation_md.append(meta_dic)
#                         data_row = []
#                         #Assign values to the different SNP types. 
#                         for indiv, variation in indiv_ids:
#                             if len(variation) == 2:
#                                 if variation == [0, 0]:
#                                     data_row.append(0)
#                                 elif variation == [0, 1]:
#                                     data_row.append(1)
#                                 elif variation == [1, 0]: 
#                                     data_row.append(1)                        
#                                 elif variation == [1, 1]:
#                                     data_row.append(2)
#                             else:
#                                 data_row.append(variation[0])
#                         data.append(data_row)
#                         previous_id = observation_id
#     if len(data) == 0:
#         raise ValueError, """No valid SNP data was present in the file. Indels will be 
# ignored"""
#     else:
#         data = array(data)
#         print 'yay'
#     return data, sample_ids, observation_ids, sample_md, observation_md
# 


#get this working
#add filter step for both bacterial and human genomic files
#human should specify ranges, 
#for instance only SNPs in all files, 

def biom_data_from_vcfs(vcfs, filter_non_snp=True, min_position=0, max_position=inf):
    oids = {}
    ordered_oids = []
    sids = {}
    ordered_sids = []
    data = {}
    for vcf in vcfs:
        for line in vcf:
            if line.startswith('##'):
                pass
            elif line.startswith('INFO'):
                pass
            elif line.startswith('#CHROM'):
                sample_ids = process_header_entry(line)
            else:
                chr_n, pos_n, rs_n, ref, alt, indiv_ids = \
                process_data_entry_line(line, sample_ids)
                try:
                    sid_index = sids[sid]
                except KeyError:
                    ordered_sids.append(sid)
                    sid_index = len(ordered_sids) - 1
                    sids[sid] = sid_index
                oid = '%s.%d' % (chrom, pos)
                if fields[4] != '.' and \
                   min_position <= pos <= max_position:
                    try:
                       oid_index = oids[oid]
                    except KeyError:
                        ordered_oids.append(oid)
                        oid_index = len(ordered_oids) - 1
                        oids[oid] = oid_index
                    # this will differ for non-haploid data:
                    data[(oid_index, sid_index)] = 1

    return data, ordered_oids, ordered_sids
    
    # data = {(0, 0): 1, (0, 1): 1, (1, 0): 1, (1, 1):0}
# oids = ['10.89673612', '10.89673554']
# #sids = [sample1, sample2]
    
def merge_otu_tables(vcf_fps):
    """Takes a list of multiple vcf files and returns a single biom table of all files."""
    master_table = None
    #open all of the files with correct extensions. Raise a value error if incorrect extension
    master_observation_ids = None
    for vcf_fp in vcf_fps:
        if vcf_fp.endswith('gz'):
            vcf_fp = gzip.open(vcf_fp)
        elif vcf_fp.endswith('vcf'):
            vcf_fp = open(vcf_fp, 'U')
        else:
            raise ValueError, "Invalid file format or extension, only '.vcf' or '.vcf.gz'\
            are accepted"
        data, sample_ids, observation_ids, sample_md, observation_md =\
        create_biom_table(vcf_fp)
        if master_observation_ids is None:
            master_observation_ids = observation_ids
        else:
            master_observation_ids = set(master_observation_ids) & set(observation_ids)
        biom_table = table_factory(data, 
                                   sample_ids, 
                                   observation_ids,
                                   sample_md, 
                                   observation_md,
                                   constructor=SparseOTUTable)
        if master_table is None:
            master_table = biom_table
        else:
            master_table.merge(biom_table)  
#         try:
#             master_table = master_table.merge(biom_table)
#         except AttributeError:
#             master_table = biom_table
    return master_table, observation_ids

def create_biom_file(vcf_fps, output_fp, mapping_fp=None, zip=None):
    master_table, master_observation_ids = merge_otu_tables(vcf_fps)
    if zip == 'gz':
        output_master_f = gzip.open(join(output_fp, 'master_table.biom.zip'), 'wb')
        output_filtered_f = gzip.open(join(output_fp, 'filtered_table.biom.zip'), 'wb')
    else:
        output_master_f = open(join(output_fp, 'master_table.biom'), 'w')
        output_filtered_f = open(join(output_fp, 'filtered_table.biom'), 'w')
    if mapping_fp is not None:
        mapping_f = MetadataMap.fromFile(mapping_fp)
        master_table.addSampleMetadata(mapping_f)
        
    master_table.getBiomFormatJsonString(generatedby(), direct_io=output_master_f)

    #create a function to filter table by
    def filter_function(values, id, md):
            return id in master_observation_ids

    filtered_table = master_table.filterObservations(filter_function)
    filtered_table.getBiomFormatJsonString(generatedby(), direct_io=output_filtered_f)
    output_master_f.close()
    output_filtered_f.close()