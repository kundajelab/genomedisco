
from __future__ import print_function
import argparse
import os
from genomedisco import concordance_utils

def main():
    command_methods = {'preprocess': concordance_utils.preprocess,
                         'concordance': concordance_utils.concordance,
                         'summary': concordance_utils.summary,
                         'cleanup':concordance_utils.clean_up,
                       'run_all': concordance_utils.run_all}
    #TODO: add command that just smooths
    command, args = concordance_utils.parse_args_genomedisco()
    if command!='cleanup':
        args['methods']='GenomeDISCO'
    command_methods[command](**args)


if __name__ == "__main__":
    main()
