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


from shutil import rmtree
from os.path import exists, join
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import remove_files, create_dir
from qiime.util import (get_qiime_temp_dir, 
                        get_tmp_filename)
from qiime.test import initiate_timeout, disable_timeout
from numpy import array


from vcf_to_biom import (process_data_entry_line, process_genotype_data, 
                         create_table_factory_objects)


class ExampleTests(TestCase):
    """ Tests the generation of lists from vcf files. """
    
    def setUp(self):
        """ Initialize variables to be used by the tests """
        self.example_file1 = example_file1 
        self.example_file2 = example_file2
        self.example_file3 = example_file3
        self.example_file4 = example_file4
        self.example_file5 = example_file5
        self.example_file6 = example_file6
#         self.example_file7 = example_file7

        
##Tests for process_data_entry_line
    def test_process_data_entry_line(self):
        """Does this work when given correct input?"""
        line = '10	89623323	rs1044322	G	A	100	PASS	.	GT:AP	0|0:0.015,0.000	0|0:0.000,0.000	0|0:0.002,0.000	0|0:0.000,0.052	0|0:0.000,0.000'
        expected = ('10', '89623323', 'rs1044322', 'G', 'A', ['0|0:0.015,0.000', '0|0:0.000,0.000', '0|0:0.002,0.000', '0|0:0.000,0.052', '0|0:0.000,0.000'])
        self.assertEqual(process_data_entry_line(line), expected)


    def test_process_data_entry_line_2(self):
        """Does this work when given correct input?"""
        line = '10	89674997	rs116819638	A	G	100	PASS	ERATE=0.0004;AN=2184;AC=13;VT=SNP;RSQ=0.9179;AA=A;AVGPOST=0.9989;LDAF=0.0063;SNPSOURCE=LOWCOV;THETA=0.0013;AF=0.01;AFR_AF=0.03	GT:DS:GL	0|0:0.000:-0.03,-1.23,-5.00	0|0:0.000:-0.02,-1.33,-5.00	0|0:0.000:-0.02,-1.33,-5.00	0|0:0.000:-0.03,-1.21,-5.00	0|0:0.000:-0.06,-0.87,-5.00	0|0:0.000:-0.48,-0.48,-0.48	0|0:0.000:-0.10,-0.69,-4.70	0|0:0.000:-0.03,-1.22,-5.00	0|0:0.000:-0.01,-1.49,-5.00	0|0:0.000:-0.00,-2.05,-5.00	0|0:0.000:-0.06,-0.87,-5.00	0|0:0.000:-0.06,-0.90,-5.00'
        expected = ('10', '89674997', 'rs116819638', 'A', 'G', ['0|0:0.000:-0.03,-1.23,-5.00', '0|0:0.000:-0.02,-1.33,-5.00', '0|0:0.000:-0.02,-1.33,-5.00', '0|0:0.000:-0.03,-1.21,-5.00', '0|0:0.000:-0.06,-0.87,-5.00', '0|0:0.000:-0.48,-0.48,-0.48', '0|0:0.000:-0.10,-0.69,-4.70', '0|0:0.000:-0.03,-1.22,-5.00', '0|0:0.000:-0.01,-1.49,-5.00', '0|0:0.000:-0.00,-2.05,-5.00', '0|0:0.000:-0.06,-0.87,-5.00', '0|0:0.000:-0.06,-0.90,-5.00'])
        self.assertEqual(process_data_entry_line(line), expected)

    def test_process_data_entry_line_3(self):
        """Does this work when given correct input?"""
        line = '10	89674997	rs116819638	A	TGC,GA	100	PASS	ERATE=0.0004;AN=2184;AC=13;VT=SNP;RSQ=0.9179;AA=A;AVGPOST=0.9989;LDAF=0.0063;SNPSOURCE=LOWCOV;THETA=0.0013;AF=0.01;AFR_AF=0.03	GT:DS:GL	0|0:0.000:-0.03,-1.23,-5.00	0|0:0.000:-0.02,-1.33,-5.00	0|0:0.000:-0.02,-1.33,-5.00	0|0:0.000:-0.03,-1.21,-5.00	0|0:0.000:-0.06,-0.87,-5.00	0|0:0.000:-0.48,-0.48,-0.48	0|0:0.000:-0.10,-0.69,-4.70	0|0:0.000:-0.03,-1.22,-5.00	0|0:0.000:-0.01,-1.49,-5.00	0|0:0.000:-0.00,-2.05,-5.00	0|0:0.000:-0.06,-0.87,-5.00	0|0:0.000:-0.06,-0.90,-5.00'
        expected = ('10', '89674997', 'rs116819638', 'A', 'TGC,GA', ['0|0:0.000:-0.03,-1.23,-5.00', '0|0:0.000:-0.02,-1.33,-5.00', '0|0:0.000:-0.02,-1.33,-5.00', '0|0:0.000:-0.03,-1.21,-5.00', '0|0:0.000:-0.06,-0.87,-5.00', '0|0:0.000:-0.48,-0.48,-0.48', '0|0:0.000:-0.10,-0.69,-4.70', '0|0:0.000:-0.03,-1.22,-5.00', '0|0:0.000:-0.01,-1.49,-5.00', '0|0:0.000:-0.00,-2.05,-5.00', '0|0:0.000:-0.06,-0.87,-5.00', '0|0:0.000:-0.06,-0.90,-5.00'])
        self.assertEqual(process_data_entry_line(line), expected)
    
    def test_process_genotype_data(self):
        """Does this work when given correct input?"""
        data = {}
        genotypes = ['0|0:0.015,0.000', '0|0:0.000,0.000', '0|0:0.002,0.000', '0|0:0.000,0.052', '0|0:0.000,0.000']
        observation_id_index = 0
        expected = {}
        self.assertEqual(process_genotype_data(data, genotypes, observation_id_index), expected)
        
    def test_process_genotype_data_1(self):
        """Does this work when given correct input?"""
        data = {}
        genotypes = ['0|0:0.015,0.000', '1|0:0.000,0.000', '0|0:0.002,0.000', '0|0:0.000,0.052', '0|0:0.000,0.000']
        observation_id_index = 1

        expected = {(1, 1):1}
        self.assertEqual(process_genotype_data(data, genotypes, observation_id_index), expected)
        
    def test_process_genotype_data_2(self):
        """Does this work when given correct input?"""
        data = {(1, 1):1}
        genotypes = ['0|0:0.000:-0.03,-1.23,-5.00', '0|0:0.000:-0.02,-1.33,-5.00', '1|1:0.000:-0.02,-1.33,-5.00', '0|0:0.000:-0.03,-1.21,-5.00', '0|1:0.000:-0.06,-0.87,-5.00', '0|0:0.000:-0.48,-0.48,-0.48', '0|0:0.000:-0.10,-0.69,-4.70', '0|0:0.000:-0.03,-1.22,-5.00', '0|0:0.000:-0.01,-1.49,-5.00', '1|0:0.000:-0.00,-2.05,-5.00', '0|0:0.000:-0.06,-0.87,-5.00', '0|0:0.000:-0.06,-0.90,-5.00']
        observation_id_index = 2
        expected = {(1, 1):1, (2, 2):2, (2,4):1, (2, 9):1}
        self.assertEqual(process_genotype_data(data, genotypes, observation_id_index), expected)

    def test_process_genotype_data_3(self):
        """Does this work when given correct input?"""
        data = {(1, 1):1, (2, 2):2, (2,4):1, (2, 9):1}
        genotypes = ['0|0:0.000:-0.03,-1.23,-5.00', '0|0:0.000:-0.02,-1.33,-5.00', '1|1:0.000:-0.02,-1.33,-5.00', '0|0:0.000:-0.03,-1.21,-5.00', '0|1:0.000:-0.06,-0.87,-5.00', '0|0:0.000:-0.48,-0.48,-0.48', '0|0:0.000:-0.10,-0.69,-4.70', '0|0:0.000:-0.03,-1.22,-5.00', '0|0:0.000:-0.01,-1.49,-5.00', '1|0:0.000:-0.00,-2.05,-5.00', '0|0:0.000:-0.06,-0.87,-5.00', '0|0:0.000:-0.06,-0.90,-5.00']
        observation_id_index = 2
        expected = {(1, 1):1, (2, 2):2, (2,4):1, (2, 9):1}
        self.assertEqual(process_genotype_data(data, genotypes, observation_id_index), expected)
        
    def test_process_genotype_data_4(self):
        """Does this work when given correct input?"""
        data = {}
        genotypes = ['0|0:0.015,0.000', '1|2:0.000,0.000', '0|0:0.002,0.000', '0|0:0.000,0.052', '0|0:0.000,0.000']
        observation_id_index = 1
        with self.assertRaises(ValueError):
            process_genotype_data(data, genotypes, observation_id_index)
            
    def test_process_genotype_data_5(self):
        """Does this work when given correct input?"""
        data = {}
        genotypes = ['0|0:0.015,0.000', '1|.:0.000,0.000', '0|0:0.002,0.000', '0|0:0.000,0.052', '0|0:0.000,0.000']
        observation_id_index = 1
        with self.assertRaises(ValueError):
            process_genotype_data(data, genotypes, observation_id_index)
                   
#Test for create_biom_table
    def test_create_table_factory_objects(self):
        """Does the function return the correct output when given the correct input?"""
        table_data = list(create_table_factory_objects(self.example_file1))
        expected = [{(0, 4):0}, ['10:89623323'], ['HG00096', 'HG00097', 'HG00099', 'HG00100', 'HG00101']]
        self.assertEqual(table_data, expected)
 
    def test_create_table_factory_objects_2(self):
        """Does the function return the correct output when given the correct input?"""
        table_data = list(create_table_factory_objects(self.example_file2))
        expected = [{(1, 11):0}, ['10:89674917', '10:89674997'], ['HG00096', 'HG00097', 'HG00099', 'HG00100', 'HG00101', 'HG00102', 'HG00103', 'HG00104', 'HG00106', 'HG00108', 'HG00109', 'HG00110']]
        self.assertEqual(table_data, expected)
        
    def test_create_table_factory_objects_3(self):
        """Does the function return the correct output when given input that is not correct?"""
        self.assertRaises(ValueError, create_table_factory_objects, self.example_file3)
        
    def test_create_table_factory_objects_4(self):
        """Does the function return the correct output when given the correct input?"""
        table_data = list(create_table_factory_objects(self.example_file4))
        expected = [{(3, 1):1},['10:89674917', '10:89674997', '10:89675036', '10:89675296'], ['HG00096', 'HG00097']]
        self.assertEqual(table_data, expected)

    def test_create_table_factory_objects_5(self):
        """Does the function return the correct output when given the correct input?"""
        table_data = list(create_table_factory_objects(self.example_file5))
        expected = [{(0, 0):2, (1, 0):1, (2, 0):1, (3, 0):1, (4, 0):1, (5, 0):1},['1:25611035', '1:25627613', '1:25627628', '1:25656531', '1:25656673', '1:25688901'], ['TNT028']]
        self.assertEqual(table_data, expected)
        
    def test_create_table_factory_objects_6(self):
        """Does the function return the correct output when given the correct input?"""
        table_data = list(create_table_factory_objects(self.example_file6))
        expected = [{(2, 6): 1, (4, 5): 1, (2, 11): 1, (2, 1): 1, (2, 13): 1, (5, 15): 0, (2, 2): 1, (2, 5): 2, (2, 4): 1}, ['10:89673416', '10:89673554', '10:89673569', '10:89673576', '10:89673612', '10:89673786'], ['HG00096', 'HG00097', 'HG00099', 'HG00100', 'HG00101', 'HG00102', 'HG00103', 'HG00104', 'HG00106', 'HG00108', 'HG00109', 'HG00110', 'HG00111', 'HG00112', 'HG00113', 'HG00114']]
        self.assertEqual(table_data, expected)
#     
example_file1 = """##fileformat=VCFv4.0													
##source=BCM:SNPTools:hapfuse													
##reference=1000Genomes-NCBI37													
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">													
##FORMAT=<ID=AP,Number=2,Type=Float,Description="Allelic Probability, P(Allele=1|Haplotype)">													
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG00096	HG00097	HG00099	HG00100	HG00101
10	89623323	rs1044322	G	A	100	PASS	.	GT:AP	0|0:0.015,0.000	0|0:0.000,0.000	0|0:0.002,0.000	0|0:0.000,0.052	0|0:0.000,0.000""".split('\n')
        
        
example_file2 = """##reference=GRCh37																				
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG00096	HG00097	HG00099	HG00100	HG00101	HG00102	HG00103	HG00104	HG00106	HG00108	HG00109	HG00110
10	89674917	rs182708158	T	G	100	PASS	ERATE=0.0004;LDAF=0.0024;AA=T;AN=2184;THETA=0.0008;VT=SNP;AVGPOST=0.9993;RSQ=0.9008;AC=5;SNPSOURCE=LOWCOV;AF=0.0023;ASN_AF=0.01	GT:DS:GL	0|0:0.000:-0.03,-1.22,-5.00	0|0:0.000:-0.01,-1.62,-5.00	0|0:0.000:-0.01,-1.87,-5.00	0|0:0.000:-0.18,-0.48,-2.19	0|0:0.000:-0.12,-0.62,-3.85	0|0:0.000:-0.14,-0.56,-3.22	0|0:0.000:-0.03,-1.15,-5.00	0|0:0.000:-0.01,-1.74,-5.00	0|0:0.000:-0.03,-1.19,-5.00	0|0:0.000:-0.05,-0.99,-5.00	0|0:0.000:-0.08,-0.75,-5.00	0|0:0.000:-0.02,-1.28,-5.00
10	89674997	rs116819638	A	G	100	PASS	ERATE=0.0004;AN=2184;AC=13;VT=SNP;RSQ=0.9179;AA=A;AVGPOST=0.9989;LDAF=0.0063;SNPSOURCE=LOWCOV;THETA=0.0013;AF=0.01;AFR_AF=0.03	GT:DS:GL	0|0:0.000:-0.03,-1.23,-5.00	0|0:0.000:-0.02,-1.33,-5.00	0|0:0.000:-0.02,-1.33,-5.00	0|0:0.000:-0.03,-1.21,-5.00	0|0:0.000:-0.06,-0.87,-5.00	0|0:0.000:-0.48,-0.48,-0.48	0|0:0.000:-0.10,-0.69,-4.70	0|0:0.000:-0.03,-1.22,-5.00	0|0:0.000:-0.01,-1.49,-5.00	0|0:0.000:-0.00,-2.05,-5.00	0|0:0.000:-0.06,-0.87,-5.00	0|0:0.000:-0.06,-0.90,-5.00""".split('\n')

example_file3 = """##fileformat=VCFv4.0													
##source=BCM:SNPTools:hapfuse													
##reference=1000Genomes-NCBI37													
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">													
##FORMAT=<ID=AP,Number=2,Type=Float,Description="Allelic Probability, P(Allele=1|Haplotype)">													
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG00096	HG00097	HG00099	HG00100	HG00101
10	89623323	rs1044322	GT	A,TC	100	PASS	.	GT:AP	1|0:0.015,0.000	1|1:0.000,0.000	0|0:0.002,0.000	0|0:0.000,0.052	1|1:0.000,0.000""".split('\n')

example_file4 = """##fileformat=VCFv4.1										
##INFO=<ID=LDAF,Number=1,Type=Float,Description=""MLE Allele Frequency Accounting for LD"">										
##INFO=<ID=AVGPOST,Number=1,Type=Float,Description=""Average posterior probability from MaCH/Thunder"">										
##INFO=<ID=RSQ,Number=1,Type=Float,Description=""Genotype imputation quality from MaCH/Thunder"">								
##INFO=<ID=ERATE,Number=1,Type=Float,Description=""Per-marker Mutation rate from MaCH/Thunder"">										
##INFO=<ID=THETA,Number=1,Type=Float,Description=""Per-marker Transition rate from MaCH/Thunder"">										
##INFO=<ID=CIEND,Number=2,Type=Integer,Description=""Confidence interval around END for imprecise variants"">										
##INFO=<ID=SNPSOURCE,Number=.,Type=String,Description=""indicates if a snp was called when analysing the low coverage or exome alignment data"">									
##reference=GRCh37										
##reference=GRCh37										
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG00096	HG00097
10	89674917	rs182708158	T	G	100	PASS	ERATE=0.0004;LDAF=0.0024;AA=T;AN=2184;THETA=0.0008;VT=SNP;AVGPOST=0.9993;RSQ=0.9008;AC=5;SNPSOURCE=LOWCOV;AF=0.0023;ASN_AF=0.01	GT:DS:GL	0|0:0.000:-0.03,-1.22,-5.00	0|0:0.000:-0.01,-1.62,-5.00
10	89674997	rs116819638	A	G	100	PASS	ERATE=0.0004;AN=2184;AC=13;VT=SNP;RSQ=0.9179;AA=A;AVGPOST=0.9989;LDAF=0.0063;SNPSOURCE=LOWCOV;THETA=0.0013;AF=0.01;AFR_AF=0.03	GT:DS:GL	0|0:0.000:-0.03,-1.23,-5.00	0|0:0.000:-0.02,-1.33,-5.00
10	89674997	rs151009112	T	C	100	PASS	ERATE=0.0005;AA=T;AN=2184;LDAF=0.0018;THETA=0.0005;VT=SNP;RSQ=0.4609;SNPSOURCE=LOWCOV;AC=2;AVGPOST=0.9974;AF=0.0009;ASN_AF=0.0035	GT:DS:GL	0|0:0.000:-0.03,-1.23,-5.00	0|0:0.000:-0.04,-1.05,-5.00
10	89675036	rs111627758	C	T	100	PASS	ERATE=0.0004;LDAF=0.0047;AA=T;AN=2184;RSQ=0.8612;THETA=0.0005;VT=SNP;AC=9;SNPSOURCE=LOWCOV;AVGPOST=0.9984;AF=0.0041;AMR_AF=0.01;EUR_AF=0.01	GT:DS:GL	0|0:0.000:-0.03,-1.16,-5.00	0|0:0.000:-0.02,-1.28,-5.00
10	89675296	rs1234224	A	G	100	PASS	ERATE=0.0004;AC=915;RSQ=0.9931;THETA=0.0004;AA=G;AN=2184;VT=SNP;LDAF=0.4188;SNPSOURCE=LOWCOV;AVGPOST=0.9961;AF=0.42;ASN_AF=0.48;AMR_AF=0.43;AFR_AF=0.42;EUR_AF=0.37	GT:DS:GL	0|0:0.000:-0.01,-1.50,-5.0	1|0:1.000:-5.00,-0.00,-4.40""".split('\n')

example_file5 = """##reference=file://E:\Genomes\Homo_sapiens\UCSC\hg19\Sequence\WholeGenomeFASTA\genome.fa									
##source=GATK 1.6									
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TNT028
1	25611035	rs2301153	G	C	10465.01	PASS	AC=2;AF=1.0;AN=2;DP=470;QD=22.27;TI=NM_001127691,NM_016124;GI=RHD,RHD;FC=Silent,Silent	GT:AD:DP:GQ:PL:VF:GQX	1/1:0,470:470:99:10465,713,0:1.000:99
1	25627613	.	C	A	4033.01	PASS	AC=1;AF=0.5;AN=2;DP=167;QD=24.15;TI=NM_001127691,NM_016124;GI=RHD,RHD;FC=Silent,Silent	GT:AD:DP:GQ:PL:VF:GQX	0/1:55,112:167:99:4033,0,891:0.671:99
1	25627627	.	C	CG	3582	PASS	AC=1;AF=0.5;AN=2;DP=163;QD=21.98;TI=NM_001127691,NM_016124;GI=RHD,RHD;FC=Noncoding,Noncoding	GT:AD:DP:GQ:PL:VF:GQX	0/1:57,106:167:99:3582,0,1837:0.650:99
1	25627628	.	A	C	3555.01	PASS	AC=1;AF=0.5;AN=2;DP=163;QD=21.81;TI=NM_001127691,NM_016124;GI=RHD,RHD;FC=Silent,Silent	GT:AD:DP:GQ:PL:VF:GQX	0/1:57,106:163:99:3555,0,982:0.650:99
1	25656531	.	T	C	2414.01	PASS	AC=1;AF=0.5;AN=2;DP=532;QD=4.54;TI=NM_001127691,NM_016124;GI=RHD,RHD;FC=Silent,Silent;EXON	GT:AD:DP:GQ:PL:VF:GQX	0/1:387,145:532:99:2444,0,7968:0.273:99
1	25656673	rs667179	C	T	2233.01	PASS	AC=1;AF=0.5;AN=2;DP=626;QD=3.57;TI=NM_001127691,NM_016124;GI=RHD,RHD;FC=Silent,Silent;EXON	GT:AD:DP:GQ:PL:VF:GQX	0/1:483,143:629:99:2263,0,9941:0.228:99
1	25688901	.	T	G	11.34	LowVariantFreq	AC=1;AF=0.5;AN=2;DP=22;QD=0.52;TI=NM_138617,NM_020485,NM_138618,NM_138616;GI=RHCE,RHCE,RHCE,RHCE;FC=Silent,Silent,Silent,Silent;EXON	GT:AD:DP:GQ:PL:VF:GQX	0/1:20,2:22:40.83:41,0,442:0.091:11""".split('\n')

example_file6 = """##INFO=<ID=ASN_AF,Number=1,Type=Float,Description="Allele Frequency for samples from ASN based on AC/AN">																								
##INFO=<ID=AFR_AF,Number=1,Type=Float,Description="Allele Frequency for samples from AFR based on AC/AN">																								
##INFO=<ID=EUR_AF,Number=1,Type=Float,Description="Allele Frequency for samples from EUR based on AC/AN">																								
##INFO=<ID=VT,Number=1,Type=String,Description="indicates what type of variant the line represents">																								
##INFO=<ID=SNPSOURCE,Number=.,Type=String,Description="indicates if a snp was called when analysing the low coverage or exome alignment data">																								
##reference=GRCh37																								
##reference=GRCh37																								
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG00096	HG00097	HG00099	HG00100	HG00101	HG00102	HG00103	HG00104	HG00106	HG00108	HG00109	HG00110	HG00111	HG00112	HG00113	HG00114
10	89673416	rs187289149	T	A	100	PASS	AA=T;AN=2184;RSQ=0.6136;LDAF=0.0018;VT=SNP;THETA=0.0006;SNPSOURCE=LOWCOV;ERATE=0.0003;AC=2;AVGPOST=0.9982;AF=0.0009;EUR_AF=0.0026	GT:DS:GL	0|0:0.000:-0.01,-1.49,-5.00	0|0:0.000:-0.01,-1.59,-5.00	0|0:0.000:-0.01,-1.54,-5.00	0|0:0.000:-0.00,-3.28,-5.00	0|0:0.000:-0.01,-1.64,-5.00	0|0:0.000:-0.19,-0.46,-2.18	0|0:0.000:-0.01,-1.49,-5.00	0|0:0.000:-0.00,-2.12,-5.00	0|0:0.000:-0.00,-2.43,-5.00	0|0:0.000:-0.00,-2.08,-5.00	0|0:0.000:-0.00,-2.74,-5.00	0|0:0.000:-0.03,-1.15,-5.00	0|0:0.000:-0.01,-1.47,-5.00	0|0:0.000:-0.05,-0.97,-5.00	0|0:0.000:-0.02,-1.39,-5.00	0|0:0.000:-0.00,-3.27,-5.00
10	89673554	rs74562518	A	C	100	PASS	ERATE=0.0004;AN=2184;VT=SNP;AVGPOST=0.9993;AA=A;SNPSOURCE=LOWCOV;AC=1;RSQ=0.6021;THETA=0.0007;LDAF=0.0008;AF=0.0005;EUR_AF=0.0013	GT:DS:GL	0|0:0.000:-0.01,-1.78,-5.00	0|0:0.000:-0.01,-1.82,-5.00	0|0:0.000:-0.01,-1.56,-5.00	0|0:0.000:-0.00,-2.70,-5.00	0|0:0.000:-0.03,-1.13,-5.00	0|0:0.000:-0.11,-0.66,-5.00	0|0:0.000:-0.00,-2.07,-5.00	0|0:0.000:-0.01,-1.81,-5.00	0|0:0.000:-0.00,-2.44,-5.00	0|0:0.000:-0.00,-2.36,-5.00	0|0:0.000:-0.03,-1.13,-5.00	0|0:0.000:-0.11,-0.65,-5.00	0|0:0.000:-0.01,-1.76,-5.00	0|0:0.000:-0.06,-0.87,-5.00	0|0:0.000:-0.11,-0.65,-4.70	0|0:0.000:-0.01,-1.54,-5.00
10	89673569	rs1234225	C	T	100	PASS	AVGPOST=0.9903;AA=T;AN=2184;VT=SNP;AC=955;THETA=0.0006;ERATE=0.0011;LDAF=0.4395;SNPSOURCE=LOWCOV;RSQ=0.9862;AF=0.44;ASN_AF=0.47;AMR_AF=0.44;AFR_AF=0.50;EUR_AF=0.37	GT:DS:GL	0|0:0.000:-0.00,-2.08,-5.00	1|0:1.000:-5.00,0.00,-5.00	1|0:1.000:-0.61,-0.12,-5.00	0|0:0.000:-0.00,-2.69,-5.00	0|1:1.000:-5.00,-0.00,-2.45	1|1:2.000:-5.00,-1.58,-0.01	0|1:1.000:-4.40,-0.00,-5.00	0|0:0.000:-0.00,-2.03,-5.00	0|0:0.000:-0.00,-2.77,-5.00	0|0:0.000:-0.00,-2.06,-5.00	0|0:0.000:-0.01,-1.64,-5.00	0|1:1.000:-0.10,-0.67,-4.70	0|0:0.000:-0.00,-2.36,-5.00	0|1:1.050:-4.70,-0.23,-0.39	0|0:0.000:-0.06,-0.88,-5.00	0|0:0.000:-0.01,-1.48,-5.00
10	89673576	rs192081057	G	A	100	PASS	AA=G;AN=2184;RSQ=0.8875;THETA=0.0005;AVGPOST=0.9996;LDAF=0.0015;VT=SNP;SNPSOURCE=LOWCOV;ERATE=0.0003;AC=3;AF=0.0014;AFR_AF=0.01	GT:DS:GL	0|0:0.000:-0.00,-2.07,-5.00	0|0:0.000:-0.00,-3.92,-5.00	0|0:0.000:-0.00,-2.57,-5.00	0|0:0.000:-0.00,-2.99,-5.00	0|0:0.000:-0.02,-1.33,-5.00	0|0:0.000:-0.01,-1.62,-5.00	0|0:0.000:-0.01,-1.77,-5.00	0|0:0.000:-0.00,-2.29,-5.00	0|0:0.000:-0.00,-2.74,-5.00	0|0:0.000:-0.01,-1.78,-5.00	0|0:0.000:-0.01,-1.66,-5.00	0|0:0.000:-0.10,-0.67,-4.70	0|0:0.000:-0.00,-2.65,-5.00	0|0:0.000:-0.11,-0.67,-4.70	0|0:0.000:-0.06,-0.88,-5.00	0|0:0.000:-0.03,-1.20,-5.00
10	89673612	rs78938505	T	C	100	PASS	RSQ=0.9920;LDAF=0.0511;AA=T;AN=2184;AC=111;THETA=0.0008;AVGPOST=0.9990;VT=SNP;SNPSOURCE=LOWCOV;ERATE=0.0003;AF=0.05;ASN_AF=0.10;AMR_AF=0.09;EUR_AF=0.03	GT:DS:GL	0|0:0.000:-0.00,-2.39,-5.00	0|0:0.000:-0.00,-2.35,-5.00	0|0:0.000:-0.00,-2.86,-5.00	0|0:0.000:-0.00,-3.05,-5.00	0|0:0.000:-0.00,-2.94,-5.00	1|0:1.000:-5.00,-0.00,-2.74	0|0:0.000:-0.01,-1.76,-5.00	0|0:0.000:-0.00,-2.02,-5.00	0|0:0.000:-0.00,-3.40,-5.00	0|0:0.000:-0.00,-2.37,-5.00	0|0:0.000:-0.01,-1.93,-5.00	0|0:0.000:-0.19,-0.47,-2.27	0|0:0.000:-0.00,-2.67,-5.00	0|0:0.000:-0.03,-1.23,-5.00	0|0:0.000:-0.02,-1.38,-5.00	0|0:0.000:-0.00,-2.62,-5.00
10	89673786	rs78068525	T	A	100	PASS	LDAF=0.0116;AA=T;AN=2184;VT=SNP;RSQ=0.9802;SNPSOURCE=LOWCOV;THETA=0.0007;ERATE=0.0003;AVGPOST=0.9995;AC=25;AF=0.01;AMR_AF=0.01;AFR_AF=0.04;EUR_AF=0.0013	GT:DS:GL	0|0:0.000:-0.05,-0.93,-5.00	0|0:0.000:-0.00,-2.34,-5.00	0|0:0.000:-0.00,-4.10,-5.00	0|0:0.000:-0.00,-2.06,-5.00	0|0:0.000:-0.01,-1.82,-5.00	0|0:0.000:-0.00,-2.85,-5.00	0|0:0.000:-0.01,-1.79,-5.00	0|0:0.000:-0.00,-4.22,-5.00	0|0:0.000:-0.00,-2.57,-5.00	0|0:0.000:-0.01,-1.78,-5.00	0|0:0.000:-0.00,-2.19,-5.00	0|0:0.000:-0.02,-1.41,-5.00	0|0:0.000:-0.00,-2.37,-5.00	0|0:0.000:-0.03,-1.24,-5.00	0|0:0.000:-0.03,-1.12,-5.00	0|0:0.000:-0.03,-1.22,-5.00""".split('\n')

# example_file6 = """##INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">									
# ##INFO=<ID=RPA,Number=.,Type=Integer,Description="Number of times tandem repeat unit is repeated, for each allele (including reference)">									
# ##INFO=<ID=RU,Number=1,Type=String,Description="Tandem repeat unit (bases)">									
# ##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">									
# ##INFO=<ID=STR,Number=0,Type=Flag,Description="Variant is a short tandem repeat">									
# ##UnifiedGenotyper="analysis_type=UnifiedGenotyper input_file=[TW10598_100_2.bam] read_buffer_size=null phone_home=STANDARD gatk_key=null tag=NA read_filter=[BadCigar] intervals=null excludeIntervals=null interval_set_rule=UNION interval_merging=ALL interval_padding=0 reference_sequence=/home/jsahl/wgfast/Ecoli/run_subsample/scratch/reference.fasta nonDeterministicRandomSeed=false disableRandomization=false maxRuntime=-1 maxRuntimeUnits=MINUTES downsampling_type=BY_SAMPLE downsample_to_fraction=null downsample_to_coverage=250 baq=OFF baqGapOpenPenalty=40.0 fix_misencoded_quality_scores=false allow_potentially_misencoded_quality_scores=false performanceLog=null useOriginalQualities=false BQSR=null quantize_quals=0 disable_indel_quals=false emit_original_quals=false preserve_qscores_less_than=6 globalQScorePrior=-1.0 allow_bqsr_on_reduced_bams_despite_repeated_warnings=false defaultBaseQualities=-1 validation_strictness=SILENT remove_program_records=false keep_program_records=false unsafe=null num_threads=2 num_cpu_threads_per_data_thread=1 num_io_threads=0 monitorThreadEfficiency=false num_bam_file_handles=null read_group_black_list=null pedigree=[] pedigreeString=[] pedigreeValidationType=STRICT allow_intervals_with_unindexed_bam=false generateShadowBCF=false logging_level=INFO log_to_file=null help=false version=false genotype_likelihoods_model=SNP pcr_error_rate=1.0E-4 computeSLOD=false annotateNDA=false pair_hmm_implementation=ORIGINAL min_base_quality_score=17 max_deletion_fraction=0.05 min_indel_count_for_genotyping=5 min_indel_fraction_per_sample=0.25 indel_heterozygosity=1.25E-4 indelGapContinuationPenalty=10 indelGapOpenPenalty=45 indelHaplotypeSize=80 indelDebug=false ignoreSNPAlleles=false allReadsSP=false ignoreLaneInfo=false reference_sample_calls=(RodBinding name= source=UNBOUND) reference_sample_name=null sample_ploidy=1 min_quality_score=1 max_quality_score=40 site_quality_prior=20 min_power_threshold_for_calling=0.95 min_reference_depth=100 exclude_filtered_reference_sites=false heterozygosity=0.001 genotyping_mode=DISCOVERY output_mode=EMIT_ALL_CONFIDENT_SITES standard_min_confidence_threshold_for_calling=30.0 standard_min_confidence_threshold_for_emitting=30.0 alleles=(RodBinding name= source=UNBOUND) max_alternate_alleles=6 contamination_fraction_to_filter=0.05 contamination_fraction_per_sample_file=null p_nonref_model=EXACT_INDEPENDENT logRemovedReadsFromContaminationFiltering=null exactcallslog=null dbsnp=(RodBinding name= source=UNBOUND) comp=[] out=org.broadinstitute.sting.gatk.io.stubs.VariantContextWriterStub no_cmdline_in_header=org.broadinstitute.sting.gatk.io.stubs.VariantContextWriterStub sites_only=org.broadinstitute.sting.gatk.io.stubs.VariantContextWriterStub bcf=org.broadinstitute.sting.gatk.io.stubs.VariantContextWriterStub debug_file=null metrics_file=null annotation=[] excludeAnnotation=[] filter_mismatching_base_and_quals=false"									
# ##contig=<ID=NC_007779.1,length=4646332>									
# ##reference=file:///home/jsahl/wgfast/Ecoli/run_subsample/scratch/reference.fasta									
# #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	vac6wt
# NC_007779.1	10098	.	C	.	68	.	AN=1;DP=1;MQ=60.00;MQ0=0	GT:DP:MLPSAC:MLPSAF	0:01
# NC_007779.1	10099	.	G	.	65	.	AN=1;DP=1;MQ=60.00;MQ0=0	GT:DP:MLPSAC:MLPSAF	0:01
# NC_007779.1	10100	.	G	.	68	.	AN=1;DP=1;MQ=60.00;MQ0=0	GT:DP:MLPSAC:MLPSAF	0:01
# NC_007779.1	10101	.	T	.	63	.	AN=1;DP=1;MQ=60.00;MQ0=0	GT:DP:MLPSAC:MLPSAF	0:01""".split('\n')
# 
# example_file7 = """##reference=file:///home/jsahl/wgfast/Ecoli/run_subsample/scratch/reference.fasta										
# #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	vac6wt	vac7wt
# INFO  13:50:05,918 SAMDataSource$SAMReaders - Initializing SAMRecords in serial 										
# INFO  13:50:05,920 SAMDataSource$SAMReaders - Done initializing BAM readers: total time 0.00 										
# INFO  13:50:06,000 AFCalcFactory - Requested ploidy 1 maxAltAlleles 6 not supported by requested model EXACT_INDEPENDENT looking for an option 										
# INFO  13:50:06,000 AFCalcFactory - Requested ploidy 1 maxAltAlleles 6 not supported by requested model EXACT_INDEPENDENT looking for an option 										
# INFO  13:50:06,000 AFCalcFactory - Selecting model EXACT_GENERAL_PLOIDY 										
# INFO  13:50:06,001 AFCalcFactory - Selecting model EXACT_GENERAL_PLOIDY 										
# NC_007779.1	10098	.	C	T	68	.	AN=1;DP=1;MQ=60.00;MQ0=0	GT:DP:MLPSAC:MLPSAF	0:01	0:01
# NC_007779.1	10099	.	G	.	65	.	AN=1;DP=1;MQ=60.00;MQ0=0	GT:DP:MLPSAC:MLPSAF	0:01	.
# NC_007779.1	10100	.	G	A	68	.	AN=1;DP=1;MQ=60.00;MQ0=0	GT:DP:MLPSAC:MLPSAF	0:01	1:01""".split('\n')
# 

if __name__ == "__main__":
    main()
    
