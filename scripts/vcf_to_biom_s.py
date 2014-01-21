#!/usr/bin/env python
# File created on 11 Jul 2013
from __future__ import division

__author__ = "John Chase"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["John Chase"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "John Chase"
__email__ = "chasejohnh@gmail.com"
__status__ = "Development"

from os.path import exists
from qiime.util import parse_command_line_parameters, make_option
from vcf_to_biom import create_biom_file
from os import makedirs


script_info = {}
script_info['brief_description'] = """Create a biom file from a vcf file"""

script_info['script_description'] = """This script takes as input a vcf file and returns
a biom formatted file. Currently some of the information in the vcf file is lost in the 
conversion"""

script_info['script_usage'] = [("Create biom file",
"Create a biom file from a vcf file",
"%prog -i input_file.vcf -o input_file.biom")]



script_info['output_description']= """The output file of this script is a biom file. In
the biom file the 'rows' are the snps identified by their locus for example 10:9892879, 
and the 'columns' are the individual ids."""

script_info['required_options'] = [\
    make_option('-i', '--input_filepaths',         
        help='Input vcf file or files. These can be a comma seperated list or contained \
in quotes with a wild card. If two or more files are passed in the output will be \
a biom table containing all of the SNPs, and a biom table containing only the \
SNPs common to all files.',
        type='existing_filepaths'),
    make_option('-o', '--output_filepath',
        help="Path to directory that will contain the output files. \
One will be created if it doesn't exist.",
        type='new_dirpath')
]
        
script_info['optional_options'] = [\
    make_option('-m', '--mapping_fp', type='existing_filepath', 
        help='Metadata mapping file filepath'),
    make_option('-z', '--zip_file', action='store_true', 
        dest='zip_file',
        help='Output will be a gzipped file'\
        ' [default: %default]'),
#     make_option('--multiple_files', action='store_true', help='Use to specify that '
#     'the -i option is composed of multiple filepaths',
#     default=False)
]
script_info['version'] = __version__

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
        
    output_fp = opts.output_filepath
    zip_file = opts.zip_file
    if zip_file:
        zip = 'gz'
        output_fp = '%s.%s' % (output_fp, zip)
    else:
        zip = None         
        
#     if exists(output_fp):
#         # don't overwrite existing output directory - make the user provide a
#         # different name or move/delete the existing directory since it may
#         # have taken a while to create.
#         option_parser.error("Output directory (%s) already exists. "
#                             "Won't overwrite." % opts.output_filepath)
        
        
        
    try:
        makedirs(output_fp)
    except OSError:
            option_parser.error("Output directory already exists. Please choose"
                " a different directory")
        
    create_biom_file(opts.input_filepaths,
                     output_fp,
                     opts.mapping_fp, 
                     zip)

if __name__ == "__main__":
    main()