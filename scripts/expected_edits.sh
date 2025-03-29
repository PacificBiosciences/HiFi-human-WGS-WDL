#!/bin/bash

# Comprehensive script to compare expected VCF with input VCF files
# This script handles both small variants and structural variants separately
# It extracts unique regions, normalizes variants, and generates a detailed report

# Set up error handling
set -e
set -o pipefail

# Define input parameters
expected_vcf="$1"
small_variants_vcf="$2"
deepvariant_vcf="$3"
ref="$4"
output_dir="$5"

# Check if required arguments are provided
if [ -z "$expected_vcf" ] || [ -z "$small_variants_vcf" ] || [ -z "$deepvariant_vcf" ] || [ -z "$output_dir" ]; then
    echo "Usage: $0    "
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$output_dir"
mkdir -p "$output_dir/small_variants_output"
mkdir -p "$output_dir/deepvariant_output"
mkdir -p "$output_dir/temp"

echo "Starting variant comparison process..."
echo "Expected VCF: $expected_vcf"
echo "Small variants VCF: $small_variants_vcf"
echo "DeepVariant VCF: $deepvariant_vcf"
echo "Output directory: $output_dir"

# Extract unique chromosome regions from expected VCF
echo "Extracting unique chromosome regions from expected VCF..."
bcftools view "$expected_vcf" | grep -v "^##" | grep -v "^#" | cut -f1 | sort | uniq > "$output_dir/temp/unique_regions.txt"
echo "Found the following unique regions:"
cat "$output_dir/temp/unique_regions.txt"

# Filter input VCFs to only include variants in those regions
echo "Filtering input VCFs to include only variants in the extracted regions..."
pattern=$(paste -sd'|' "${output_dir}/temp/unique_regions.txt")
grep -E '^#|^('"$pattern"')' "$small_variants_vcf" > "$output_dir/temp/filtered_small_variants.vcf"
grep -E '^#|^('"$pattern"')' "$deepvariant_vcf" > "$output_dir/temp/filtered_deepvariant.vcf"


# bcftools view -R "$output_dir/temp/unique_regions.txt" "$small_variants_vcf" > "$output_dir/temp/filtered_small_variants.vcf"
# bcftools view -R "$output_dir/temp/unique_regions.txt" "$deepvariant_vcf" > "$output_dir/temp/filtered_deepvariant.vcf"

# Normalize VCFs
echo "Normalizing VCFs..."




bcftools norm  -m -any -f ${ref} "$expected_vcf" -O v | bcftools sort -O v -o "${output_dir}/temp/normalized_expected.vcf"
bgzip "${output_dir}/temp/normalized_expected.vcf"
bcftools index "${output_dir}/temp/normalized_expected.vcf.gz"

bcftools norm  -m -any -f ${ref} "${output_dir}/temp/filtered_small_variants.vcf" -O v | bcftools sort -O v -o "${output_dir}/temp/normalized_small_variants.vcf"
bgzip "${output_dir}/temp/normalized_small_variants.vcf"
bcftools index "${output_dir}/temp/normalized_small_variants.vcf.gz"

bcftools norm  -m -any -f ${ref} "${output_dir}/temp/filtered_deepvariant.vcf" -O v | bcftools sort -O v -o "${output_dir}/temp/normalized_deepvariant.vcf"
bgzip "${output_dir}/temp/normalized_deepvariant.vcf"
bcftools index "${output_dir}/temp/normalized_deepvariant.vcf.gz"


# Run bcftools isec to find intersections
echo "Running bcftools isec to find intersections..."
bcftools isec -p "$output_dir/small_variants_output" "$output_dir/temp/normalized_expected.vcf.gz" "$output_dir/temp/normalized_small_variants.vcf.gz"
bcftools isec -p "$output_dir/deepvariant_output" "$output_dir/temp/normalized_expected.vcf.gz" "$output_dir/temp/normalized_deepvariant.vcf.gz"

echo "Intersection analysis complete. Generating report..."

export output_dir

python << 'EOF'
import pandas as pd 
import os 
import sys
import glob

def read_vcf(file_path):
    """Read a VCF file into a pandas DataFrame."""
    try:
        # Skip header lines that start with '##'
        with open(file_path, 'r') as f:
            header_lines = 0
            for line in f:
                if line.startswith('##'):
                    header_lines += 1
                else:
                    break
        
        # Read the file, skipping the header lines
        df = pd.read_csv(
            file_path,
            sep='\t',
            skiprows=header_lines,
        )
        
        # Check if DataFrame is empty
        if df.empty:
            #print("The resulting DataFrame is empty.")
            return pd.DataFrame()
            
        
        # Ensure the column names are correct by removing leading '#' from columns
        if df.columns[0].startswith('#'):
            df.columns = [col.lstrip('#') for col in df.columns]
        
        return df
    
    except Exception as e:
        print(f"Error reading VCF file {file_path}: {e}")
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
    report = []

    # Process expected-only variants
    if os.path.exists(expected_only_path):
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
                            "Status": "Expected only - not found in input",
                            "Allele_Presence": "N/A"
                        })
        except Exception as e:
             print(f"Error processing expected-only variants for {variant_type}: {e}")

    # Process shared variants
    if os.path.exists(shared_path):
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
                    "Status": "Found in both files - MATCH",
                    "Allele_Presence": allele_status
                })
                print(report)
        except Exception as e:
            print(f"Error processing shared variants for {variant_type}: {e}")

    return report


#Main execution
output_dir = os.environ['output_dir']
print(f"Using output directory: {output_dir}")

#Process both variant types
small_variants_report = process_intersection_results(output_dir, "small_variants")
deepvariant_report = process_intersection_results(output_dir, "deepvariant")

#Combine reports
combined_report = small_variants_report + deepvariant_report

#Create report DataFrame and save to CS
if combined_report:
    report_df = pd.DataFrame(combined_report)
    report_path = os.path.join(output_dir, "variant_comparison_report.csv")
    report_df.to_csv(report_path, index=False)
    print(f"Report generated: {report_path}")

    # Generate summary statistics
    summary = {
        "Total variants in expected VCF": len(small_variants_report) + len(deepvariant_report),
        "Small variants found in both files": len([v for v in small_variants_report if v["Status"] == "Found in both files - MATCH"]),
        "Small variants only in expected VCF": len([v for v in small_variants_report if v["Status"] == "Expected only - not found in input"]),
        "DeepVariant variants found in both files": len([v for v in deepvariant_report if v["Status"] == "Found in both files - MATCH"]),
        "DeepVariant variants only in expected VCF": len([v for v in deepvariant_report if v["Status"] == "Expected only - not found in input"])
    }

    # Print summary
    print("\nSummary Statistics:")
    for key, value in summary.items():
        print(f"{key}: {value}")
else:
    print('No variants found for comparison')
EOF


echo "Variant comparison process completed."
echo "Results are available in $output_dir/variant_comparison_report.csv"

# Optional: Clean up temporary files
# rm -rf "$output_dir/temp"


