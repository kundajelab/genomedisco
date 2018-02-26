
from __future__ import print_function
import argparse
import os
from genomedisco import genomedisco_utils

def main():
    command_methods = {'split': genomedisco_utils.split_by_chromosome,
                         'reproducibility': genomedisco_utils.reproducibility,
                         'summary': genomedisco_utils.summary,
                       'run_all': genomedisco_utils.run_all}
    command, args = genomedisco_utils.parse_args()
    global repo_dir
    global bashrc_file

    repo_dir=os.path.dirname(os.path.realpath(__file__))
    bashrc_file=repo_dir+'/configuration_files/bashrc.configuration'

    command_methods[command](**args)


if __name__ == "__main__":
    main()
