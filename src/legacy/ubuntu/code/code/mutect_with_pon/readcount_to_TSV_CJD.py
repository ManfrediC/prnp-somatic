#!/usr/bin/env python3
import sys
import csv
import os
import glob

def parse_bam_readcount(input_file, output_file):
    """
    Convert bam-readcount output to a readable TSV format
    """
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t', lineterminator='\n')
        
        header = [
            "CHROM", "POS", "REF", "BASE", "COUNT", "MEAN_BQ", "MEAN_MQ", 
            "MEAN_POS", "FWD", "REV", "FRAC", "MMLQ", "MEAN_READ_POS", 
            "FWD_BQ", "REV_BQ", "FRD", "SB_RATIO"
        ]
        writer.writerow(header)
        
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
                
            chrom, pos, ref, depth = row[0], row[1], row[2], row[3]
            
            for allele_field in row[4:]:
                parts = allele_field.split(':')
                base = parts[0]
                metrics = parts[1:]
                
                if len(metrics) < 13:
                    metrics += [''] * (13 - len(metrics))
                
                out = [
                    chrom, pos, ref, base,
                    metrics[0],  # COUNT
                    metrics[1],  # MEAN_BQ
                    metrics[2],  # MEAN_MQ
                    metrics[3],  # MEAN_POS
                    metrics[4],  # FWD
                    metrics[5],  # REV
                    metrics[6],  # FRAC
                    metrics[7],  # MMLQ
                    metrics[8],  # MEAN_READ_POS
                    metrics[9],  # FWD_BQ
                    metrics[10], # REV_BQ
                    metrics[11], # FRD
                    metrics[12], # SB_RATIO
                ]
                writer.writerow(out)

def process_directory(input_dir, output_dir):
    """
    Process all .txt files in the input directory and save results to output directory
    """
    os.makedirs(output_dir, exist_ok=True)
    input_files = glob.glob(os.path.join(input_dir, "*.txt"))
    
    for input_file in input_files:
        base_filename = os.path.basename(input_file)
        sample_name = os.path.splitext(base_filename)[0]
        output_filename = f"{sample_name}_metrics.tsv"
        output_path = os.path.join(output_dir, output_filename)
        
        parse_bam_readcount(input_file, output_path)
        print(f"Processed {input_file} -> {output_path}")

def main():
    """
    Main function to handle processing of sample files only
    """
    sample_input_dir = "/add/results/with_PoN/readcounts/"
    sample_output_dir = "/add/results/with_PoN/metrics/"
    
    print("Processing samples (CJD and controls)...")
    process_directory(sample_input_dir, sample_output_dir)
    
    print("All files processed successfully!")

if __name__ == "__main__":
    main()
