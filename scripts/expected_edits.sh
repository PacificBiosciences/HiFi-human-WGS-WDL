#!/bin/bash

# Comprehensive script to compare expected VCF with input VCF files
# This script handles both small variants and structural variants separately
# It extracts unique regions, normalizes variants, and generates a detailed report

# Set up error handling
set -e -o pipefail

# Define input parameters
expected_vcf="$1"
small_variants_vcf="$2"
structural_variants_vcf="$3"
ref="$4"
output_dir="$5"

# Check if required arguments are provided
if [ -z "$expected_vcf" ] || [ -z "$small_variants_vcf" ] || [ -z "$structural_variants_vcf" ] || [ -z "$output_dir" ]; then
    echo "Usage: $0 <expected.vcf> <small.vcf> <structural_variants.vcf> <ref.fa> <output_dir>"
    exit 1
fi
# Create output directory if it doesn't exist
mkdir -p "$output_dir"
mkdir -p "$output_dir/small_variants_output"
mkdir -p "$output_dir/structural_variants_output"
mkdir -p "$output_dir/temp"




echo "Starting variant comparison process..." >&2
echo "Expected VCF: $expected_vcf" >&2
echo "Small variants VCF: $small_variants_vcf" >&2
echo "Structural variants VCF: $structural_variants_vcf" >&2
echo "Output directory: $output_dir" >&2

#function to check if vcf file is bgzipped
#if it is, then original file name is returned
#if not, the file is then comprezzed and the new file name (with .gz appended) is returned
check_bgzip() {
    local vcf_file="$1"
    if [[ "$vcf_file" == *.gz ]]; then
        echo "$vcf_file"
    else
        echo "File $vcf_file is not bgzipped. Compressing ..." >&2
        bgzip -c "$vcf_file" > ${output_dir}/temp/${vcf_file}.gz
        bcftools index ${output_dir}/temp/${vcf_file}.gz
        echo "${output_dir}/temp/${vcf_file}.gz"
    fi
}

#checking if each vcf file input is zipped and setting new variables to the bgzipped versions
expected_vcf_bgz=$(check_bgzip "$expected_vcf")
small_variants_vcf_bgz=$(check_bgzip "$small_variants_vcf")
structural_variants_vcf_bgz=$(check_bgzip "$structural_variants_vcf")

echo "Using the following bgzipped VCFs:" >&2
echo "Expected VCF: $expected_vcf_bgz" >&2
echo "Small Variants VCF: $small_variants_vcf_bgz" >&2
echo "Structural Variants VCF: $structural_variants_vcf_bgz" >&2

echo "Sorting expected edits" >&2
bcftools sort "${expected_vcf_bgz}" -O b -o "${output_dir}/temp/expected_edits_sorted.vcf.gz" 

echo "Indexing expected edits" >&2
bcftools index "${output_dir}/temp/expected_edits_sorted.vcf.gz"

# Filter input VCFs to only include variants in those regions
echo "Filtering input VCFs to include only variants in the extracted regions..." >&2
bcftools view -R "${output_dir}/temp/expected_edits_sorted.vcf.gz" --regions-overlap pos "$small_variants_vcf_bgz" -O b > "$output_dir/temp/filtered_small_variants.vcf.gz"
bcftools view -R "${output_dir}/temp/expected_edits_sorted.vcf.gz" --regions-overlap pos "$structural_variants_vcf_bgz" -O b > "$output_dir/temp/filtered_structural_variants.vcf.gz"

# Normalize VCFs
echo "Normalizing VCFs..." >&2

bcftools norm  -m -any -f ${ref} "$expected_vcf_bgz" -O b | bcftools sort -O b -o "${output_dir}/temp/normalized_expected.vcf.gz"
bcftools index "${output_dir}/temp/normalized_expected.vcf.gz"

bcftools norm  -m -any -f ${ref} "${output_dir}/temp/filtered_small_variants.vcf.gz" -O b | bcftools sort -O b -o "${output_dir}/temp/normalized_small_variants.vcf.gz"
bcftools index "${output_dir}/temp/normalized_small_variants.vcf.gz"

bcftools norm  -m -any -f ${ref} "${output_dir}/temp/filtered_structural_variants.vcf.gz" -O b | bcftools sort -O b -o "${output_dir}/temp/normalized_structural_variants.vcf.gz"
bcftools index "${output_dir}/temp/normalized_structural_variants.vcf.gz"


# Run bcftools isec to find intersections
echo "Running bcftools isec to find intersections..." >&2
bcftools isec -p "$output_dir/small_variants_output" "$output_dir/temp/normalized_expected.vcf.gz" "$output_dir/temp/normalized_small_variants.vcf.gz"
bcftools isec -p "$output_dir/structural_variants_output" "$output_dir/temp/normalized_expected.vcf.gz" "$output_dir/temp/normalized_structural_variants.vcf.gz"

echo "Intersection analysis complete. Generating report..." >&2

export output_dir

python << 'EOF'
import pandas as pd 
import os 
import sys
import glob

def read_vcf(file_path):
    """Read a VCF file into a pandas DataFrame."""
    try:
        with open(file_path, 'r') as f:
            header_count = 0
            for header_count, line in enumerate(f):
                if line.startswith('#CHROM'):
                    break
            cols = line.lstrip("#").strip().split("\t")
        df=pd.read_csv(file_path, sep="\t", skiprows=header_count+1, header=None, names=cols)        
        print(df)
        return df
    
    except Exception as e:
        print(f"Error reading VCF file {file_path}: {e}", file=sys.stderr)
        return pd.DataFrame()

def check_allele_presence(gt_field):
    """Check if variant is present on both, one, or no alleles"""
    if not isinstance(gt_field, str):
        return "Unknown"
    gt = gt_field.split(':') if ':' in gt_field else gt_field

    if gt in ["1/1", "1|1"]:
        return "Both alleles"
    elif gt in ["0/1", "1/0", "0|1", "1|0"]:
        return "One allele"
    elif gt in ["0/0", "0|0"]:
        return "No alleles"
    else:
        return "Complex/Other"

def process_intersection_results(output_dir, variant_type):
    """Process the intersection results for a specific variant type"""
    # Define file paths
    expected_only_path = os.path.join(output_dir, f"{variant_type}_output", "0000.vcf")
    input_only_path = os.path.join(output_dir, f"{variant_type}_output", "0001.vcf")
    shared_path = os.path.join(output_dir, f"{variant_type}_output", "0002.vcf")

    # Check that all required files exist
    required_files = {
        "Expected-only VCF": expected_only_path,
        "Input-only VCF": input_only_path,
        "Shared VCF": shared_path
    }
    for desc, file_path in required_files.items():
        if not os.path.exists(file_path):
            print(f"Error: Required file '{desc}' does not exist: {file_path}. Aborting processing.", file=sys.stderr)
            sys.exit(1)

    report = []

    # Process expected-only variants
    try:
        expected_only = read_vcf(expected_only_path)
        if not expected_only.empty: 
            for _, row in expected_only.iterrows():
                    report.append({
                        "Variant_Type": variant_type,
                        "CHROM": row["CHROM"],
                        "POS": row["POS"],
                        "REF": row["REF"],
                        "ALT": row["ALT"],
                        "Status": "NOT CONFIRMED in sample",
                        "Allele_Presence": "N/A"
                    })
    except Exception as e:
         print(f"Error processing expected-only variants for {variant_type}: {e}", file=sys.stderr)

    # Process shared variants
    try:
        shared = read_vcf(shared_path)
        #print(shared.columns)
        for _, row in shared.iterrows():
            # Get genotype field if available (assuming it's in the FORMAT and sample columns)
            #print(row)
            gt_field = "Unknown"
            if len(row) > 9:  # If there are sample columns
                format_col = row["FORMAT"] if "FORMAT" in row else row.iloc
                #print(format_col)
                sample_col = row.iloc[9] # First sample column
                #print(row.iloc[9])
                #print(sample_col)
                if isinstance(format_col, str) and isinstance(sample_col, str):
                    format_fields = format_col.split(':')
                    sample_fields = sample_col.split(':')
                    
                    if 'GT' in format_fields:
                        gt_index = format_fields.index('GT')
                        if gt_index < len(sample_fields):
                            gt_field = sample_fields[gt_index]
            #print(gt_field)
            allele_status = check_allele_presence(gt_field)
            
            report.append({
                "Variant_Type": variant_type,
                "CHROM": row["CHROM"],
                "POS": row["POS"],
                "REF": row["REF"],
                "ALT": row["ALT"],
                "Status": "CONFIRMED in sample",
                "Allele_Presence": allele_status
            })
            print(report)
    except Exception as e:
        print(f"Error processing shared variants for {variant_type}: {e}", file=sys.stderr)

    return report


#Main execution
output_dir = os.environ['output_dir']
print(f"Using output directory: {output_dir}")

#Process both variant types
small_variants_report = process_intersection_results(output_dir, "small_variants")
structural_variants_report = process_intersection_results(output_dir, "structural_variants")

#Combine reports
combined_report = small_variants_report + structural_variants_report

#Create report DataFrame and save to CS
if combined_report:
    report_df = pd.DataFrame(combined_report)
    report_path = os.path.join(output_dir, "variant_comparison_report.csv")
    report_df.to_csv(report_path, index=False)
    print(f"Report generated: {report_path}")

    # Generate summary statistics
    summary = {
        "Total edit confirmations": len(small_variants_report) + len(structural_variants_report),
        "Edits confirmed by small variants": len([v for v in small_variants_report if v["Status"] == "CONFIRMED in sample"]),
        "Edits NOT confirmed by small variants": len([v for v in small_variants_report if v["Status"] == "NOT CONFIRMED in sample"]),
        "Edits confirmed by structural variants": len([v for v in structural_variants_report if v["Status"] == "CONFIRMED in sample - MATCH"]),
        "Edits NOT confirmed by structural variants": len([v for v in structural_variants_report if v["Status"] == "NOT CONFIRMED in sample"])
    }

    # Print summary
    print("\nSummary Statistics:")
    for key, value in summary.items():
        print(f"{key}: {value}")
else:
    print('No variants found for comparison', file=sys.stderr)
EOF


# # Pass the output directory to the Python script
# python3 - "$output_dir"
#
echo "Variant comparison process completed." >&2
echo "$output_dir/variant_comparison_report.csv"

# Optional: Clean up temporary files
# rm -rf "$output_dir/temp"


