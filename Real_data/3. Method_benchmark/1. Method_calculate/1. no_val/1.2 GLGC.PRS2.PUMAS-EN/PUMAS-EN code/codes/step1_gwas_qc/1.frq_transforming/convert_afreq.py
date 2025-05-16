#!/usr/bin/env python
"""
Convert 1000 Genomes afreq format to standard frequency format.

This script converts frequency files from the format used in 1000 Genomes to
a standard format with columns: CHR, SNP, A1, A2, MAF, NCHROBS.
The script ensures MAF is always ≤ 0.5 by swapping alleles if needed.
"""

import pandas as pd
import argparse
import os

def convert_afreq_to_standard(input_file, output_file):
    """
    Convert a .afreq file to standard frequency format with MAF ≤ 0.5.
    
    Parameters:
    -----------
    input_file : str
        Path to input .afreq file
    output_file : str
        Path to save the converted output file
    """
    print(f"Converting {input_file} to standard format...")
    
    # Load the afreq file
    df = pd.read_csv(input_file, delim_whitespace=True)
    
    # Rename columns to match standard format
    df.rename(columns={
        '#CHROM': 'CHR', 
        'ID': 'SNP', 
        'REF': 'A1', 
        'ALT': 'A2', 
        'ALT_FREQS': 'MAF', 
        'OBS_CT': 'NCHROBS'
    }, inplace=True)
    
    # Make a copy of original A1/A2
    df['A1_orig'] = df['A1'].copy()
    df['A2_orig'] = df['A2'].copy()
    
    # Identify SNPs where ALT_FREQS >= 0.5
    mask = df['MAF'] >= 0.5
    
    # For SNPs where ALT_FREQS >= 0.5, swap REF and ALT and adjust MAF
    df.loc[mask, 'A1'] = df.loc[mask, 'A2_orig']
    df.loc[mask, 'A2'] = df.loc[mask, 'A1_orig']
    df.loc[mask, 'MAF'] = 1 - df.loc[mask, 'MAF']
    
    # For SNPs where ALT_FREQS < 0.5, no swap needed but we need to reorder columns
    # From REF=A1, ALT=A2 to REF=A2, ALT=A1
    not_mask = ~mask
    df.loc[not_mask, 'A1'] = df.loc[not_mask, 'A2_orig']
    df.loc[not_mask, 'A2'] = df.loc[not_mask, 'A1_orig']
    
    # Drop temporary columns
    df.drop(['A1_orig', 'A2_orig'], axis=1, inplace=True)
    
    # Round MAF to 4 decimal places
    df['MAF'] = df['MAF'].round(4)
    
    # Select and reorder columns to match the standard format
    df = df[['CHR', 'SNP', 'A1', 'A2', 'MAF', 'NCHROBS']]
    
    # Save the output file
    df.to_csv(output_file, sep='\t', index=False)
    
    print(f'Converted file saved to {output_file}')
    print(f'Number of SNPs processed: {len(df)}')
    print(f'Number of SNPs with MAF adjustment: {mask.sum()}')

def main():
    parser = argparse.ArgumentParser(description='Convert afreq files to standard frequency format')
    parser.add_argument('--input_file', required=True, help='Path to input .afreq file')
    parser.add_argument('--output_file', required=True, help='Path to save the converted output file')
    
    args = parser.parse_args()
    
    output_dir = os.path.dirname(args.output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    convert_afreq_to_standard(args.input_file, args.output_file)

if __name__ == "__main__":
    main()