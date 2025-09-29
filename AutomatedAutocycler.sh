#!/bin/bash

# This script runs the full AutoCycler pipeline on all FASTQ and BAM files
# in the reads/ directory. It organises the output into step-specific main
# directories, with each containing subdirectories for every isolate.

# --- Configurable Parameters ---
FALLBACK_GENOME_SIZE="5000000" # Fallback genome size if auto-detection fails
THREADS=120 # Number of CPU threads to use. CHANGE THIS to match your system.
READ_TYPE=pacbio_hifi # Read type [default: ont_r10] [possible values: ont_r9, ont_r10, pacbio_clr, pacbio_hifi]
READS_DIR="reads" # Path to directory containing input read files
checkm2_database=/home/mcknid01/softwareDependencies/CheckM2_database/uniref100.KO.1.dmnd  # CheckM2 database path

# --- Initialize conda for this script ---
# This allows conda activate to work within the script
eval "$(conda shell.bash hook)"
conda activate autocycler

# --- Create Top-Level Step Directories Before the Loop ---
echo "Creating main output directories."
mkdir -p 00_filtered_reads 01_subsampled_reads 02_assemblies 03_autocycler_out 04_racon_polish 05_final_consensus

# --- Main Loop ---
# Find all common FASTQ and BAM files and loop through them.
for reads_file in "${READS_DIR}"/*.fastq.gz "${READS_DIR}"/*.fastq "${READS_DIR}"/*.fq.gz "${READS_DIR}"/*.fq "${READS_DIR}"/*.bam; do
    # If no files are found, the loop will run once with the literal string.
    # This check skips that iteration if the file doesn't actually exist.
    [ -f "$reads_file" ] || continue

    # 1. Set up variables for this specific file
    # Remove extensions to get the base name
    base_name=$(basename "${reads_file}" .fastq.gz | sed 's/\.fastq$//' | sed 's/\.fq\.gz$//' | sed 's/\.fq$//' | sed 's/\.bam$//')

    echo "================================================================="
    echo "Starting AutoCycler pipeline for: ${reads_file}"
    echo "================================================================="

    # Define the output directories for this specific isolate within each step
    step0_dir="00_filtered_reads/${base_name}"
    step1_dir="01_subsampled_reads/${base_name}"
    step2_dir="02_assemblies/${base_name}"
    step3_dir="03_autocycler_out/${base_name}"
    step4_dir="04_racon_polish/${base_name}"

    # Convert BAM to fastq if needed
    if [[ "${reads_file}" == *.bam ]]; then
        # Define the path for the converted FASTQ file
        mkdir -p "${step0_dir}"
        converted_fastq="${step0_dir}/${base_name}_converted.fastq"

        # Check if the converted FASTQ file already exists
        if [ -f "${converted_fastq}" ]; then
            echo "Detected existing converted FASTQ file. Skipping BAM conversion."
            reads_file_for_filtering="${converted_fastq}"
        else
            echo "Detected BAM file. Converting to FASTQ using samtools."
            samtools fastq -@ "${THREADS}" "${reads_file}" > "${converted_fastq}"
            reads_file_for_filtering="${converted_fastq}"
        fi
    else
        reads_file_for_filtering="${reads_file}"
    fi

    # Step 0a: Filter reads with filtlong
    echo "Step 0a: Filtering reads with filtlong into ${step0_dir}."
    mkdir -p "${step0_dir}"
    filtered_fastq="${step0_dir}/${base_name}_trimmed.fastq"

    if [ -f "${filtered_fastq}" ]; then
        echo " 	 	Existing filtered reads detected. Skipping filtlong."
        reads_file_for_autocycler="${filtered_fastq}"
    else
        echo " 	 	Running filtlong (min_length=1000, keep_percent=95)..."
        filtlong --min_length 1000 --keep_percent 95 "${reads_file_for_filtering}" > "${filtered_fastq}"
        reads_file_for_autocycler="${filtered_fastq}"
    fi

    # Step 0b: Automatically determine genome size
    echo "Step 0b: Determining genome size from filtered reads."
    genome_size_file="${step0_dir}/${base_name}_genome_size.txt"

    if [ -f "${genome_size_file}" ]; then
        echo " 	 	Existing genome size file detected. Reading cached value."
        genome_size=$(cat "${genome_size_file}")
        echo " 	 	Using cached genome size: ${genome_size}"
    else
        echo " 	 	Running genome size estimation..."
        genome_size=$(autocycler helper genome_size --reads "${reads_file_for_autocycler}" --threads "${THREADS}")
        
        if [ -z "${genome_size}" ]; then
            echo " 	 	Warning: Could not automatically determine genome size. Using fallback value: ${FALLBACK_GENOME_SIZE}"
            genome_size="${FALLBACK_GENOME_SIZE}"
        else
            echo " 	 	Detected genome size: ${genome_size}"
            # Cache the genome size for future runs
            echo "${genome_size}" > "${genome_size_file}"
        fi
    fi

    # Step 1: Subsampling
    echo "Step 1: Subsampling reads into ${step1_dir}."
    if [ -d "${step1_dir}" ] && compgen -G "${step1_dir}/sample_*.fastq" > /dev/null; then
        echo " 	 	Existing subsamples detected. Skipping subsampling"
    else
        autocycler subsample --reads "${reads_file_for_autocycler}" --out_dir "${step1_dir}" --genome_size "${genome_size}"
    fi

    # Step 2: Assembling subsamples
    echo "Step 2: Assembling subsamples into ${step2_dir}."
    if [ -d "${step2_dir}" ] && [ "$(ls -A "${step2_dir}" 2>/dev/null)" ]; then
        echo " 	 	Existing assemblies detected. Skipping assembly."
    else
        mkdir -p "${step2_dir}"
        for assembler in canu flye metamdbg miniasm necat nextdenovo plassembler raven; do
            for i in 01 02 03 04; do
                autocycler helper "${assembler}" \
                    --reads "${step1_dir}/sample_${i}.fastq" \
                    --out_prefix "${step2_dir}/${assembler}_${i}" \
                    --read_type "${READ_TYPE}" \
                    --threads "${THREADS}" \
                    --genome_size "${genome_size}"
            done
        done
    fi
  
    # Step 3: Compressing assemblies
    echo "Step 3: Compressing assemblies into ${step3_dir}."
    if [ -d "${step3_dir}" ] && [ "$(ls -A "${step3_dir}" 2>/dev/null)" ]; then
        echo " 	 	Existing compressed output detected. Skipping compression of contigs."
    else
        autocycler compress -i "${step2_dir}" -a "${step3_dir}"
    fi

    # Step 4: Clustering contigs
    echo "Step 4: Clustering contigs."
    cluster_dir="${step3_dir}/clustering"
    if [ -d "${cluster_dir}" ] && [ "$(ls -A "${cluster_dir}" 2>/dev/null)" ]; then
        echo " 	 	Existing clustering detected. Skipping contig clustering."
    else
        autocycler cluster -a "${step3_dir}"
    fi

    # Steps 5 & 6: Trimming and resolving clusters
    echo "Steps 5 & 6: Trimming and resolving clusters."
    qc_pass_dir="${step3_dir}/clustering/qc_pass"
    if [ -d "${qc_pass_dir}" ] && [ "$(ls -A "${qc_pass_dir}" 2>/dev/null)" ]; then
        for c in "${qc_pass_dir}"/cluster_*; do
            [ -e "$c" ] || continue
            [ -d "$c" ] || continue
            echo " 	 	 	- Processing cluster $(basename "$c")."

            # Check for existing trim results - using ls to check for any files/dirs starting with 3_trim
            if ls "$c"/2_trim* >/dev/null 2>&1; then
                echo " 	 	 	Existing trim results detected. Skipping trim."
            else
                autocycler trim -c "$c"
            fi

            # Check for existing resolve results - check for the 5_final.gfa file that resolve produces
            if [ -f "$c/5_final.gfa" ]; then
                echo " 	 	 	Existing resolve results detected. Skipping resolve."
            else
                autocycler resolve -c "$c"
            fi
        done

        # Step 7: Combine
        echo "Step 7: Combining resolved clusters into final assembly."
        combined_file="${step3_dir}/consensus_assembly.fasta"
        if [ -f "${combined_file}" ]; then
            echo " 	 	Existing combined output detected. Skipping combining."
        else
            echo " 	 	Combining resolved clusters..."
            autocycler combine -a "${step3_dir}" -i ${qc_pass_dir}/cluster_*/5_final.gfa
        fi
    fi

    # Step 8: Racon polish
    echo "Step 8: Polishing consensus with Racon."
    consensus_assembly="${step3_dir}/consensus_assembly.fasta"
    racon_output="${step4_dir}/genome.racon.fasta"

    if [ -f "${consensus_assembly}" ]; then
        if [ -f "${racon_output}" ]; then
            echo " 	 	Existing Racon polished assembly detected. Skipping polishing."
        else
            echo " 	 	Running Racon on consensus assembly."
            
            # Create the Racon working directory
            mkdir -p "${step4_dir}"
            
            # Determine minimap2 options based on read type
            if [[ "${READ_TYPE}" == "pacbio_hifi" ]] || [[ "${READ_TYPE}" == "pacbio_clr" ]]; then
                minimap2_opts="-ax map-pb"
            else
                minimap2_opts="-ax map-ont"
            fi
            
            # Align reads to consensus assembly
            echo " 	 	Aligning reads to consensus assembly with minimap2..."
            alignment_file="${step4_dir}/reads_to_assembly.sam"
            minimap2 ${minimap2_opts} -t "${THREADS}" "${consensus_assembly}" "${filtered_fastq}" > "${alignment_file}"
            
            # Run Racon
            echo " 	 	Running Racon polishing..."
            racon -t "${THREADS}" "${filtered_fastq}" "${alignment_file}" "${consensus_assembly}" > "${racon_output}"
            
            # Clean up alignment file to save space
            rm -f "${alignment_file}"
        fi
    else
        echo "Warning: Consensus assembly not found at ${consensus_assembly}. Skipping Racon."
    fi

    # Step 9: Copy and rename the polished consensus assembly to 05_final_consensus
    echo "Step 9: Copying final polished assembly to 05_final_consensus."
    final_assembly_dest="05_final_consensus/${base_name}.fasta"

    if [ -f "${racon_output}" ]; then
        if [ -f "${final_assembly_dest}" ]; then
            echo " 	 	Final assembly already exists in 05_final_consensus. Skipping copy."
        else
            echo " 	 	Copying ${racon_output} to ${final_assembly_dest}"
            cp "${racon_output}" "${final_assembly_dest}"
        fi
    elif [ -f "${consensus_assembly}" ]; then
        echo "Warning: Racon output not found. Copying unpolished consensus as fallback."
        if [ ! -f "${final_assembly_dest}" ]; then
            cp "${consensus_assembly}" "${final_assembly_dest}"
        fi
    else
        echo "Warning: No assembly found to copy to final destination."
    fi
done

# Step 10: Run CheckM2 for quality assessment
echo "Step 10: Assessing assembly quality with CheckM2."
checkm2_input_dir="05_final_consensus"
checkm2_output_dir="05_final_consensus/checkM2report"

# Check if CheckM2 has already been run for this isolate by looking for its report file
if [ -f "${checkm2_output_dir}/quality_report.tsv" ]; then
    echo "Existing CheckM2 output found. Skipping quality assessment."
else
    # Activate the CheckM2 Conda environment
    echo "Activating 'checkm2' conda environment."
    conda activate checkm2
    echo "Running CheckM2..."

    # Run CheckM2
    checkm2 predict \
        --database_path $checkm2_database \
        --input "${checkm2_input_dir}" \
        --output-directory "${checkm2_output_dir}" \
        --threads "${THREADS}" \
        -x fasta \
        --force 
        
        
    # Deactivate the environment after use
    conda deactivate
fi

echo "--- Finished pipeline for ${reads_file} ---"
echo

echo "================================================================="
echo "All files have been processed."
echo "================================================================="
