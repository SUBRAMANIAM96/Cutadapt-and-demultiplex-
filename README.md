# Cutadapt-and-demultiplex-
script for cleaning multiple primer pairs for a singke sample and tag based demultiplication using Cutadapt
#!/bin/bash

# Define the file paths
file_path="/mnt/c/Users/User/Downloads/Master thesis"
forward_read="${file_path}/Monk seal_diet pool3 F . fastq"
reverse_read="${file_path}/Monk seal_diet pool3 R . fastq"

# Define quality filtering options
quality_cutoff=20  # Minimum quality score for trimming
min_length=50      # Minimum read length after trimming

# Define primer sequences and tags for each sample (2 sets of primers per sample)
declare -A forward_primers_1
declare -A reverse_primers_1
declare -A forward_primers_2
declare -A reverse_primers_2
declare -A sample_tags

# Add the primer sequences and tags for each sample (two primer pairs per sample)
forward_primers_1["Seal_dietDes08a"]="GTCGGTAAAACTCGTGCCAGC"
reverse_primers_1["Seal_dietDes08a"]="CATAGTGGGGTATCTAATCCCAGTTTG"
forward_primers_2["Seal_dietDes08a"]="GTTGGTAAATCTCGTGCCAGC"
reverse_primers_2["Seal_dietDes08a"]="CATAGTGGGGTATCTAATCCTAGTTTG"
sample_tags["Seal_dietDes08a"]="ggtaag"

forward_primers_1["Seal_dietDes10a"]="GTCGGTAAAACTCGTGCCAGC"
reverse_primers_1["Seal_dietDes10a"]="CATAGTGGGGTATCTAATCCCAGTTTG"
forward_primers_2["Seal_dietDes10a"]="GTTGGTAAATCTCGTGCCAGC"
reverse_primers_2["Seal_dietDes10a"]="CATAGTGGGGTATCTAATCCTAGTTTG"
sample_tags["Seal_dietDes10a"]="cactct"

forward_primers_1["Seal_dietDes11a"]="GTCGGTAAAACTCGTGCCAGC"
reverse_primers_1["Seal_dietDes11a"]="CATAGTGGGGTATCTAATCCCAGTTTG"
forward_primers_2["Seal_dietDes11a"]="GTTGGTAAATCTCGTGCCAGC"
reverse_primers_2["Seal_dietDes11a"]="CATAGTGGGGTATCTAATCCTAGTTTG"
sample_tags["Seal_dietDes11a"]="aacgcg"

forward_primers_1["Seal_dietDes12a"]="GTCGGTAAAACTCGTGCCAGC"
reverse_primers_1["Seal_dietDes12a"]="CATAGTGGGGTATCTAATCCCAGTTTG"
forward_primers_2["Seal_dietDes12a"]="GTTGGTAAATCTCGTGCCAGC"
reverse_primers_2["Seal_dietDes12a"]="CATAGTGGGGTATCTAATCCTAGTTTG"
sample_tags["Seal_dietDes12a"]="aacaac"

forward_primers_1["Seal_dietDes13b1"]="GTCGGTAAAACTCGTGCCAGC"
reverse_primers_1["Seal_dietDes13b1"]="CATAGTGGGGTATCTAATCCCAGTTTG"
forward_primers_2["Seal_dietDes13b1"]="GTTGGTAAATCTCGTGCCAGC"
reverse_primers_2["Seal_dietDes13b1"]="CATAGTGGGGTATCTAATCCTAGTTTG"
sample_tags["Seal_dietDes13b1"]="taacat"

# Check if Cutadapt is installed
if ! command -v cutadapt &> /dev/null; then
    echo "Error: Cutadapt is not installed. Please install it using:"
    echo "  pip install cutadapt"
    exit 1
fi

# Check if input files exist
if [[ ! -f "$forward_read" || ! -f "$reverse_read" ]]; then
    echo "Error: One or more required files are missing!"
    [[ ! -f "$forward_read" ]] && echo "Missing: $forward_read"
    [[ ! -f "$reverse_read" ]] && echo "Missing: $reverse_read"
    exit 1
fi

# Loop through the samples and process each one
for sample in "${!forward_primers_1[@]}"; do
    forward_primer_1="${forward_primers_1[$sample]}"
    reverse_primer_1="${reverse_primers_1[$sample]}"
    forward_primer_2="${forward_primers_2[$sample]}"
    reverse_primer_2="${reverse_primers_2[$sample]}"
    tag="${sample_tags[$sample]}"

    echo "Processing sample: $sample"

    # Define output files for trimmed reads
    trimmed_forward="${file_path}/trimmed_${sample}_F.fastq"
    trimmed_reverse="${file_path}/trimmed_${sample}_R.fastq"
    
    # Step 1: Trim using both primer pairs (forward and reverse for two sets of primers)
    cutadapt -q $quality_cutoff --minimum-length $min_length \
             --pair-filter=any \
             -g "${forward_primer_1}" -G "${reverse_primer_1}" \
             -a "${forward_primer_1}" -A "${reverse_primer_1}" \
             -g "${forward_primer_2}" -G "${reverse_primer_2}" \
             -a "${forward_primer_2}" -A "${reverse_primer_2}" \
             -o "$trimmed_forward" -p "$trimmed_reverse" \
             "$forward_read" "$reverse_read"

    # Check if cutadapt ran successfully
    if [[ $? -ne 0 ]]; then
        echo "Error: Trimming failed for $sample."
        continue
    fi

    # Step 2: Remove any remaining adapter sequences
    cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC \
             -o "${file_path}/no_adapters_${sample}_F.fastq" \
             -p "${file_path}/no_adapters_${sample}_R.fastq" \
             "$trimmed_forward" "$trimmed_reverse"

    # Step 3: Demultiplex based on the tag using cutadapt for paired-end reads
    cutadapt -g "$tag" -G "$tag" \
             -o "${file_path}/demux_${sample}_F.fastq" \
             -p "${file_path}/demux_${sample}_R.fastq" \
             "${file_path}/no_adapters_${sample}_F.fastq" "${file_path}/no_adapters_${sample}_R.fastq"

    # Check if demultiplexing was successful
    if [[ $? -eq 0 ]]; then
        echo "Sample $sample processed successfully."
    else
        echo "Error processing sample $sample with demultiplexing."
    fi
done

echo "Trimming, adapter removal, and demultiplexing completed."
