import os
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description='Merge files from input directories.')
parser.add_argument('--target_dir', type=str, default='./inputs/', help='Target directory for merged files')
parser.add_argument('--input_dirs', nargs='+', required=True, help='List of input directories')
parser.add_argument('--suffix', nargs='+', required=True, help='List of suffixes for each input directory')
parser.add_argument('--file_to_merge', nargs='+', default=['AnalysisResults'], help='List of files to merge')

args = parser.parse_args()

# Extract arguments
target_dir = args.target_dir
input_dirs = args.input_dirs
suffix = args.suffix
file_to_merge = args.file_to_merge

# Loop over the input directories and suffixes
for i, (input_dir, suf) in enumerate(zip(input_dirs, suffix)):
    train_number = input_dir.split('/')[-2]
    for file in file_to_merge:
        os.system(f'alien.py cp -T 64 alien://{input_dir}/{file}.root file:{target_dir}/{file}_{suf}.root')
