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
from make_snp_significance_histogram import (get_observation_significances, get_bins,
get_observation_counts, get_observation_ratios)


class ExampleTests(TestCase):
    """ Tests the generation of lists from vcf files. """
    
    def setUp(self):
        """ Initialize variables to be used by the tests """
        self.example_file1 = example_file1
        self.example_file2 = example_file2
        self.example_file3 = example_file3


        
##Tests for process_data_entry_line
    def test_get_observation_significances(self):
        """Does this work when given correct input?"""
        sig_file = list(self.example_file1)
        expected = {89673569:2.18334244957, 89673786:1.49353654761, 
                    89673554:0.925925925958, 89673554:0.925925925958, 
                    89673612:0.372504190703, 89673576:0.261714337513, 
                    89673416:0.00296568141117}
        self.assertEqual(get_observation_significances(sig_file, 0, 1), expected)
        
    def test_get_observation_significances_2(self):
        """Does this work when given correct input?"""
        sig_file = list(self.example_file2)
        expected = {37727329: 0.000410679807173, 17122884: 1.65149355352e-05, 
                    45665159: 0.000458525064284, 20657577: 1.36374430037e-06, 
                    25487245: 0.000531817951012, 20657552: 3.67213902751e-08, 
                    34837524: 0.000278191479631, 17866198: 0.000277108761292, 
                    21681115: 6.27910706602e-10, 25636117: 0.000434635239138}
        self.assertEqual(get_observation_significances(sig_file, 0, 1), expected)
        
    def test_get_observation_significances_3(self):
        """Does this work when given correct input?"""
        sig_file = list(self.example_file3)
        expected = {89673569:2.18334244957, 89673786:0}
        self.assertEqual(get_observation_significances(sig_file, 0, 1), expected)
    
    def test_get_bins(self):
        input = {37727329: 0.000410679807173, 17122884: 1.65149355352e-05, 
                    45665159: 0.000458525064284, 20657577: 1.36374430037e-06, 
                    25487245: 0.000531817951012, 20657552: 3.67213902751e-08, 
                    34837524: 0.000278191479631, 17866198: 0.000277108761292, 
                    21681115: 6.27910706602e-10, 25636117: 0.000434635239138}
        bin_size = 5000000
        expected = [17122884, 22122884, 27122884, 32122884, 37122884, 42122884, 45665159]
        self.assertEqual(get_bins(input, bin_size, False), expected) 
    
    def test_get_bins_2(self):
        input = {89673569:2.18334244957, 89673786:1.49353654761, 
                 89673554:0.925925925958, 89673554:0.925925925958, 
                 89673612:0.372504190703, 89673576:0.261714337513, 
                 89673416:0.00296568141117}
        bin_size = 100
        expected = [89673416, 89673516, 89673616, 89673716, 89673786]
        self.assertEqual(get_bins(input, bin_size, False), expected)
        
    def test_get_bins_3(self):
        input = {89673569:2.18334244957, 89673786:1.49353654761, 
                 89673554:0.925925925958, 89673554:0.925925925958, 
                 89673612:0.372504190703, 89673576:0.261714337513, 
                 89673416:0.00296568141117}
        bin_size = 1000
        expected = [89673416, 89673786]
        self.assertEqual(get_bins(input, bin_size), expected)
        
    def test_get_bins_3(self):
        input = {89673569:2.18334244957, 89673786:1.49353654761, 
                 89673554:0.925925925958, 89673554:0.925925925958, 
                 89673612:0.372504190703, 89673576:0.261714337513, 
                 89673416:0.00296568141117}
        bin_size = 1000
        expected = [89673416, 89673786]
        self.assertEqual(get_bins(input, bin_size, False), expected)
        
        
    def test_get_bins_4(self):
        input = {89673569:2.18334244957, 89673786:1.49353654761, 
                 89673554:0.925925925958, 89673612:0.372504190703, 
                 89673576:0.261714337513, 89673416:0.00296568141117}
        bin_size = 2
        expected = [89673416, 89673569, 89673612, 89673786]
#         (89673416, 89673554, 89673569, 89673576, 89673612, 89673786]
        
        self.assertEqual(get_bins(input, bin_size, True), expected)
    
    def test_get_observation_counts(self):
        input = {37727329: 0.000410679807173, 17122884: 1.65149355352e-05, 
                 45665159: 0.000458525064284, 20657577: 1.36374430037e-06, 
                 25487245: 0.000531817951012, 20657552: 3.67213902751e-08, 
                 34837524: 0.000278191479631, 17866198: 0.000277108761292, 
                 21681115: 6.27910706602e-10, 25636117: 0.000434635239138}
        bins = [17122884, 22122884, 27122884, 32122884, 37122884, 42122884, 45665159]
        alpha = .05
        expected = ({37122884: 1, 22122884: 2, 17122884: 5, 27122884: 0, 42122884: 0, 32122884: 1, 45665159: 1}, 
                    {37122884: 1, 22122884: 2, 17122884: 5, 27122884: 0, 42122884: 0, 32122884: 1, 45665159: 1})
        self.assertEqual(get_observation_counts(input, bins, alpha), expected)
        
        
    def test_get_observation_counts_2(self):
        input = {37727329: 1, 17122884: 1.65149355352e-05, 
                 45665159: .05, 20657577: .9, 
                 25487245: 0.000531817951012, 20657552: 3.67213902751e-08, 
                 34837524: 0.000278191479631, 17866198: 0.000277108761292, 
                 21681115: 6.27910706602e-10, 25636117: 0.000434635239138}
        bins = [17122884, 22122884, 27122884, 32122884, 37122884, 42122884, 45665159]
        alpha = .05
        expected = ({17122884: 5, 22122884: 2, 27122884: 0, 32122884: 1, 37122884: 1, 42122884: 0, 45665159: 1}, 
                    {17122884: 4, 22122884: 2, 27122884: 0, 32122884: 1, 37122884: 0, 42122884: 0, 45665159: 1})
        self.assertEqual(get_observation_counts(input, bins, alpha), expected)

  #   def test_get_observation_ratios(self):
#         file = list(self.example_file1)
#         bin_size = 100
#         expected = {89673416: 20.0, 89673616: 1, 89673786: 0.0, 89673516: 0.0, 89673716: 1}
#         self.assertEqual(get_observation_ratios(file, bin_size), expected)
#         
# 
#     def test_get_observation_ratios_2(self):
#         file = list(self.example_file2)
#         bin_size = 5000000
#         expected = ({27122884: 1, 45665159: 20.0, 17122884: 20.0, 42122884: 1, 37122884: 20.0, 22122884: 20.0, 32122884: 20.0}, [17122884, 22122884, 27122884, 32122884, 37122884, 42122884, 45665159])
#         self.assertEqual(get_observation_ratios(file, bin_size), expected)
#             
        
example_file1 = """OTU	Test-Statistic	P	FDR_P	Bonferroni_P	male_mean	female_mean
10:89673569	2.18334244957	0.139511186868	0.665005662746	0.837067121208	0.841904761905	0.904761904762
10:89673786	1.49353654761	0.221668554249	0.665005662746	1.0	0.0171428571429	0.0282186948854
10:89673554	0.925925925958	0.335923813143	0.671847626287	1.0	0.0	0.00176366843034
10:89673612	0.372504190703	0.541642423898	0.730733784978	0	0.106666666667	0.0970017636684
10:89673576	0.261714337513	0.608944820815	0.730733784978	.01	0.00190476190476	0.00352733686067
10:89673416	0.00296568141117	0.956570214112	0.956570214112	1.0	0.00190476190476	0.00176366843034""".split('\n')

example_file2 = """OTU	prob	Bonferroni_corrected	FDR_corrected	male_mean	female_mean
22:21681115	6.27910706602e-10	3.96814450144e-05	3.96814450144e-05	7.34604866163e-06	4.18173940692e-06
22:20657552	3.67213902751e-08	0.00232064497982	0.00116032248991	1.12390397164e-05	7.78798320245e-06
22:20657577	1.36374430037e-06	0.0861831848061	0.0287277282687	1.0041379102e-05	7.18332367965e-06
22:17122884	1.65149355352e-05	1.04367786608	0.26091946652	1.90543735743e-05	2.18128725006e-05
22:17866198	0.000277108761292	17.5121652786	2.9186942131	2.29583227962e-05	2.07435635739e-05
22:34837524	0.000278191479631	17.5805887467	2.51151267811	8.83255860521e-06	6.66029442493e-06
22:37727329	0.000410679807173	25.9533210941	3.24416513676	1.35514318331e-05	1.59799530432e-05
22:25636117	0.000434635239138	27.4672085726	3.05191206362	1.13446880327e-05	1.37529591141e-05
22:45665159	0.000458525064284	28.9769499625	2.89769499625	5.32042101336e-06	3.64386180703e-06
22:25487245	0.000531817951012	33.6087672321	3.05534247565	6.54708719021e-06	8.59133387292e-06""".split('\n')
        
example_file3 = """OTU	prob	Bonferroni_corrected	FDR_corrected	male_mean	female_mean
10:89673569	2.18334244957	0.139511186868	0.665005662746	0.837067121208	0.841904761905	0.904761904762
10:89673786	0	0.221668554249	0.665005662746	1.0	0.0171428571429	0.0282186948854""".split('\n')

if __name__ == "__main__":
    main()