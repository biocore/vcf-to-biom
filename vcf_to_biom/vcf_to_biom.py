#!/usr/bin/env python
# File created on 11 Jul 2013
from __future__ import division

__author__ = "John Chase"
__copyright__ = "Copyright 2013, The QIIME project"
__credits__ = ["John Chase"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "John Chase"
__email__ = "chasejohnh@gmail.com"
__status__ = "Development"

from biom.table import table_factory, SparseOTUTable
from biom.parse import MetadataMap, generatedby
from numpy import array
import os
import gzip

def process_data_entry_line(line):
    '''Takes a data line from a vcf file and returns the different data types as a list'''
    fields = line.split('\t') 
    chr_n = fields[0] 
    pos_n = fields[1]
    rs_n = fields[2]
    ref = fields[3] 
    alt = fields[4]
    #genotypes still contains all of the data about the SNP call including quality 
    #information
    genotypes = fields[9:]
    return (chr_n, pos_n, rs_n, ref, alt, genotypes)
    
def process_genotype_data(data, genotypes, observation_id_index):
    '''Takes a dictionary, empty or not, and returns a dictionary. The key the coordinates
    of observation and sample id and the value is a 1, 2, or 3'''
    sample_id_index = 0
    for genotype in genotypes:
        try:
            #This parses the data so that only the genotype information is left
            genotype = map(int,genotype.split(':')[0].split('|'))
        except ValueError:
            genotype = map(int,genotype.split(':')[0].split('/'))
        if genotype == [0, 0]:
            pass
        elif genotype == [0, 1]:
            data[(observation_id_index, sample_id_index)] = 1
        elif genotype == [1, 0]:
            data[(observation_id_index, sample_id_index)] = 1
        elif genotype == [1, 1]:
            data[(observation_id_index, sample_id_index)] = 2
        else:
            raise ValueError("Can't handle genotype %s (sample id: %s, observation id: %s)" % 
                             (genotype, sample_id_index, observation_id_index))
        sample_id_index += 1
    return data

def create_table_factory_objects(vcf_file):
    ordered_observations_ids = []
    observations_ids = {}
    master_observation_ids = set([])
    data = {}
    for line in vcf_file:
        if line.startswith('##'):
            pass
        elif line.startswith('INFO'):
            pass
        elif line.startswith('#CHROM'): 
            ordered_sample_ids = line.strip().split('\t')[9:] 
        else:
            chr_n, pos_n, rs_n, ref, alt, genotypes = process_data_entry_line(line)
            observation_id = chr_n + ':' + pos_n
            if observation_id in master_observation_ids or len(ref) > 1 or len(alt) > 1:
                pass
            else:
                master_observation_ids.add(observation_id)
                ordered_observations_ids.append(observation_id)
                observation_id_index = len(ordered_observations_ids) - 1
                data = process_genotype_data(data, genotypes, observation_id_index)
    if len(ordered_observations_ids) == 0:
        raise ValueError, "No valid SNP data was found, check VCF file."
    try:
        data[(observation_id_index, (len(ordered_sample_ids) - 1))]
    except KeyError:
        data[(observation_id_index, (len(ordered_sample_ids) - 1))] = 0
    return data, ordered_observations_ids, ordered_sample_ids

# data = {(0, 0): 1, (0, 1): 1, (1, 0): 1, (1, 1):0}
# oids = ['10.89673612', '10.89673554']
# #sids = [sample1, sample2]

def create_biom_file(vcf_fp, output_fp, mapping_fp=None, zip=None):
    if vcf_fp.endswith('gz'):
        vcf_f = gzip.open(vcf_fp)
    elif vcf_fp.endswith('vcf'):
        vcf_f = open(vcf_fp, 'U')
    else:
        raise ValueError, "Invalid file format or extension, only '.vcf' or '.vcf.gz' are\
accepted"
    data, observation_ids, sample_ids =\
    create_table_factory_objects(vcf_f)
    sample_md = None
    table = table_factory(data, sample_ids=sample_ids, observation_ids=observation_ids, constructor=SparseOTUTable)
    if mapping_fp != None:
        mapping_f = MetadataMap.fromFile(mapping_fp)
        biom_table.addSampleMetadata(mapping_f)
    if zip == 'gz':
        output_f = gzip.open('%s.%s' % (output_fp, zip), 'wb')
    else:
        output_f = open(output_fp, 'w')
    table.getBiomFormatJsonString(generatedby(), direct_io=output_f)
    output_f.close()
    

                       
