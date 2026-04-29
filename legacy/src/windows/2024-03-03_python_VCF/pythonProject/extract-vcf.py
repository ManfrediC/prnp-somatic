import gzip
import sys
import re






# Function to parse INFO field and extract AD
def parse_info(info_str):
    info_dict = {}
    info_fields = info_str.split(";")
    for field in info_fields:
        key_value = field.split("=")
        if len(key_value) == 2:
            info_dict[key_value[0]] = key_value[1]
    return info_dict.get("AD", ".")

# Function to parse FORMAT and sample data and extract GT and AD
def parse_format(format_str, sample_data):
    format_fields = format_str.split(":")
    gt_index = format_fields.index("GT")
    ad_index = format_fields.index("AD")
    sample_fields = sample_data.split(":")
    gt = sample_fields[gt_index]
    ad = sample_fields[ad_index] if ad_index < len(sample_fields) else "."
    return gt, ad

# Regular expression pattern to extract sample name from VCF header
sample_pattern = re.compile(r'FORMAT\t(.+)$')

# Input VCF file path
vcf_file_path = "/add/2024-03-03_count_path_variants/CJD30_E200K_variants.vcf.gz"

# Output TSV file path
output_tsv_file = "output.tsv"

# Open gzipped VCF file
with gzip.open(vcf_file_path, "rt") as vcf_file, open(output_tsv_file, "w") as tsv_file:
    # Write header
    tsv_file.write("Sample\tCHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tGT\tPL\tAD_REF\tAD_ALT\n")

    sample_name = None

    # Parse VCF file line by line
    for line in vcf_file:
        if line.startswith("##bcftools_callCommand"):
            # Extract sample name from VCF header
            sample_match = sample_pattern.search(line)
            if sample_match:
                sample_name = sample_match.group(1)

        elif not line.startswith("#"):
            fields = line.strip().split("\t")
            chrom = fields[0]
            pos = fields[1]
            vid = fields[2]
            ref = fields[3]
            alt = fields[4]
            qual = fields[5]
            filt = fields[6]
            info = fields[7]
            format_str = fields[8]
            sample_data = fields[9]

            # Parse format and sample data
            gt, ad = parse_format(format_str, sample_data)

            # Write data to TSV file
            tsv_file.write(f"{sample_name}\t{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t{qual}\t{filt}\t{gt}\t.\t{ad}\n")

print(f"Output written to {output_tsv_file}")