#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

function log {
    local timestamp=$(date +"%Y-%m-%d %T")
    echo "---> [ ${timestamp} ] $@" >&2
}

function check_file {
    local file=$1
    if [ ! -f ${file} ]; then
        echo "ERROR: File ${file} does not exist."
        exit 1
    fi
    if [ ! -r ${file} ]; then
        echo "ERROR: File ${file} is not readable."
        exit 1
    fi
    if [ ! -s ${file} ]; then
        echo "ERROR: File ${file} is empty."
        exit 1
    fi
}

usage() { 
    echo "Usage: $0 -l <Library_Name> -i <Input_File_TSV> -o <Output_Directory> -c <Threads>" 1>&2
    exit 1
}

# Preset Values #
threads=1
out_dir=$(pwd)

while getopts ":l:i:o:c:h" opt; do
    case "${opt}" in
        l) library_name=${OPTARG} ;;
        i)
            input_tsv=${OPTARG}
            if [ ! -f ${input_tsv} ]; then echo "ERROR: The input TSV file does not exist."; usage; fi 
            ;;
        o) out_dir=${OPTARG} ;;
        c)
            threads=${OPTARG}
            if ! [[ ${threads} =~ ^[0-9]+$ ]]; then echo "ERROR: Threads must be a positive integer"; usage; fi
            ;;
        h)
            usage
            ;;
        :)
            printf "ERROR: Missing argument for -%s\n" "$OPTARG"
            usage
            ;;
        \?)
            printf "ERROR: Illegal option: -%s\n" "${OPTARG}"
            usage
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

# Variable Checks
if [ -z ${library_name+x} ]; then echo "ERROR: Missing Library Name."; usage; fi

# Parse through input TSV file
check_file ${input_tsv}
# Checking if the library_name exists in the input TSV
if ! grep -q ${library_name} ${input_tsv}; then
    echo "ERROR: Library name ${library_name} not found in the input TSV file."
    exit 1
fi

log "Starting Lindsey Analysis for ${library_name}."

working_dir=$(dirname $(realpath ${input_tsv}))
log "ALL RAW DATA Inputs should be placed in: ${working_dir}/${library_name}/"
# The input files should be in the same directory as the input TSV in their library

# Extracting the reference and fastq files from the input TSV
reference=${working_dir}/${library_name}/$(grep ${library_name} ${input_tsv} | cut -f 2)
check_file ${reference}
initial_fastq1=${working_dir}/${library_name}/$(grep ${library_name} ${input_tsv} | cut -f 3)
check_file ${fastq1}
initial_fastq2=${working_dir}/${library_name}/$(grep ${library_name} ${input_tsv} | cut -f 4)
check_file ${fastq2}
final_fastq1=${working_dir}/${library_name}/$(grep ${library_name} ${input_tsv} | cut -f 5)
check_file ${fastq1}
final_fastq2=${working_dir}/${library_name}/$(grep ${library_name} ${input_tsv} | cut -f 6)
check_file ${fastq2}

# Extracting other parameters from the input TSV
fwd_adapter=$(grep ${library_name} ${input_tsv} | cut -f 7 | tr -d '[:space:]')
fwd_rc_adapter=$(grep ${library_name} ${input_tsv} | cut -f 8 | tr -d '[:space:]')
rev_adapter=$(grep ${library_name} ${input_tsv} | cut -f 9 | tr -d '[:space:]')
rev_rc_adapter=$(grep ${library_name} ${input_tsv} | cut -f 10 | tr -d '[:space:]')

adapter1="${fwd_adapter}...${rev_rc_adapter}"
adapter2="${rev_adapter}...${fwd_rc_adapter}"

syn_codon=$(grep ${library_name} ${input_tsv} | cut -f 11 | tr -d '[:space:]')
wt_codon=$(grep ${library_name} ${input_tsv} | cut -f 12 | tr -d '[:space:]')
library_dna_sequence=$(grep ${library_name} ${input_tsv} | cut -f 13 | tr -d '[:space:]')
library_aa_sequence=$(grep ${library_name} ${input_tsv} | cut -f 14 | tr -d '[:space:]')

log "Input TSV: ${input_tsv}"
log "Reference: ${reference}"
log "Number of Cores: ${threads}"
log "Initial FASTQ_R1: ${initial_fastq1}"
log "Initial FASTQ_R2: ${initial_fastq2}"
log "Final FASTQ_R1: ${final_fastq1}"
log "Final FASTQ_R2: ${final_fastq2}"
log "Output Directory: ${out_dir}"
log "Script Directory: ${SCRIPT_DIR}"
log "Current Directory: $(pwd)"
log "Starting analysis at $(date)"
log "Creating output directory if it does not exist."
mkdir -p ${out_dir}

log "Creating temporary directory."
tmp_dir="${out_dir}/${library_name}_$(date +%s)"
mkdir -p "${tmp_dir}"
if [ ! -d ${tmp_dir} ]; then echo "ERROR: Failed to create temporary directory."; exit 1; fi

log "Building Bowtie2 index."
bowtie2-build ${reference} ${tmp_dir}/ref_index
if [ $? -ne 0 ]; then echo "ERROR: Failed to build Bowtie2 index."; exit 1; fi
log "Bowtie2 index built successfully."

# Part 1 - Alignment
log "Running Bowtie2 alignment."
bowtie2 -p ${threads} -x ${tmp_dir}/ref_index -1 ${initial_fastq1} -2 ${initial_fastq2} -S ${tmp_dir}/${library_name}_initial_aligned.sam
bowtie2 -p ${threads} -x ${tmp_dir}/ref_index -1 ${final_fastq1} -2 ${final_fastq2} -S ${tmp_dir}/${library_name}_final_aligned.sam
if [ $? -ne 0 ]; then echo "ERROR: Bowtie2 alignment failed."; exit 1; fi
log "Bowtie2 alignment completed successfully."

# Part 2 - BAM Quality Check
log "Converting SAM to BAM."
samtools view -@ ${threads} -Sb -f2 -q42 ${tmp_dir}/${library_name}_initial_aligned.sam > ${out_dir}/${library_name}_initial_aligned_mapped_MAPQ42.bam
samtools view -@ ${threads} ${tmp_dir}/${library_name}_initial_aligned.sam | cut -f 5 | sort | uniq -c | sort -n | awk '{printf("MAPQ:%s\t%d\n",$2,$1);}'
samtools view -@ ${threads} -Sb -f2 -q42 ${tmp_dir}/${library_name}_final_aligned.sam > ${out_dir}/${library_name}_final_aligned_mapped_MAPQ42.bam
samtools view -@ ${threads} ${tmp_dir}/${library_name}_final_aligned.sam | cut -f 5 | sort | uniq -c | sort -n | awk '{printf("MAPQ:%s\t%d\n",$2,$1);}'
if [ $? -ne 0 ]; then echo "ERROR: BAM quality check failed."; exit 1; fi
log "BAM quality check completed successfully."

# Part 3 - Converting BAM to FASTQ
log "Converting BAM to FASTQ."
# I changed this from bedtools bamtofastq to samtools fastq because bedtools does not support multithreading
# Samtools finished in ~1 Minute
# Bedtools finished in ~10 Minutes
samtools fastq -@ ${threads} -1 ${tmp_dir}/${library_name}_initial_mapped_MAPQ42_R1.fastq -2 ${tmp_dir}/${library_name}_initial_mapped_MAPQ42_R2.fastq -N ${out_dir}/${library_name}_initial_aligned_mapped_MAPQ42.bam
samtools fastq -@ ${threads} -1 ${tmp_dir}/${library_name}_final_mapped_MAPQ42_R1.fastq -2 ${tmp_dir}/${library_name}_final_mapped_MAPQ42_R2.fastq -N ${out_dir}/${library_name}_final_aligned_mapped_MAPQ42.bam
if [ $? -ne 0 ]; then echo "ERROR: BAM to FASTQ conversion failed."; exit 1; fi
log "BAM to FASTQ conversion completed successfully."

# Part 4 - pandaseq to merge
log "Running Pandaseq to merge FASTQ files."
pandaseq -B -f ${tmp_dir}/${library_name}_initial_mapped_MAPQ42_R1.fastq -r ${tmp_dir}/${library_name}_initial_mapped_MAPQ42_R2.fastq -w ${out_dir}/${library_name}_initial_mapped_MAPQ42_merged.fasta -g ${tmp_dir}/pandaseq.log -T ${threads}
pandaseq -B -f ${tmp_dir}/${library_name}_final_mapped_MAPQ42_R1.fastq -r ${tmp_dir}/${library_name}_final_mapped_MAPQ42_R2.fastq -w ${out_dir}/${library_name}_final_mapped_MAPQ42_merged.fasta -g ${tmp_dir}/pandaseq.log -T ${threads}
if [ $? -ne 0 ]; then echo "ERROR: Pandaseq merging failed."; exit 1; fi
log "Pandaseq merging completed successfully."

# Part 5 - Trimming
log "Running cutadapt for trimming."
cutadapt -j ${threads} -a ${adapter1} -a ${adapter2} -o ${out_dir}/${library_name}_initial_mapped_MAPQ42_merged_trimmed.fasta ${out_dir}/${library_name}_initial_mapped_MAPQ42_merged.fasta
cutadapt -j ${threads} -a ${adapter1} -a ${adapter2} -o ${out_dir}/${library_name}_final_mapped_MAPQ42_merged_trimmed.fasta ${out_dir}/${library_name}_final_mapped_MAPQ42_merged.fasta
if [ $? -ne 0 ]; then echo "ERROR: Cutadapt trimming failed."; exit 1; fi
log "Cutadapt trimming completed successfully."

# Part 6 - Python Code
log "Running Python script for further analysis."
python3 ${SCRIPT_DIR}/uniq_ALC1_codons.py -s ${out_dir}/${library_name}_initial_mapped_MAPQ42_merged_trimmed -l ${library_name} --sequence ${library_dna_sequence} --syn ${syn_codon} --wt ${wt_codon}
python3 ${SCRIPT_DIR}/uniq_ALC1_codons.py -s ${out_dir}/${library_name}_final_mapped_MAPQ42_merged_trimmed -l ${library_name} --sequence ${library_dna_sequence} --syn ${syn_codon} --wt ${wt_codon}
python3 ${SCRIPT_DIR}/logfitness_ALC1.py -s1 ${out_dir}/${library_name}_final_mapped_MAPQ42_merged_trimmed.WT_spike.RCPM -s2 ${out_dir}/${library_name}_initial_mapped_MAPQ42_merged_trimmed.WT_spike.RCPM -l ${library_name} --sequence ${library_aa_sequence} --output ${out_dir}
if [ $? -ne 0 ]; then echo "ERROR: Python script execution failed."; exit 1; fi

# Check to see if the output files were created
if [ ! -f ${out_dir}/diflog_${library_name}_final_mapped_MAPQ42_merged_trimmed.WT_spike.RCPM_${library_name}_initial_mapped_MAPQ42_merged_trimmed.WT_spike.RCPM.csv ]; then
    echo "ERROR: Output file not created."
    exit 1
fi
log "Python script executed successfully."

# Remove temporary files
log "Cleaning up temporary files."
rm -rf ${tmp_dir}
if [ $? -ne 0 ]; then echo "ERROR: Failed to remove temporary files."; exit 1; fi
log "Temporary files removed successfully."
log "Analysis completed successfully."
log "Output files are located in: ${out_dir}"
log "Script completed at $(date)"
log "Script finished successfully."
log "All done!"
