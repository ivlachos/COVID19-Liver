
TISSUE=$1
CORES=$2
FASTQPATH=$3



# cutadapt settings
trimQ=10 #quality to trim 3'
minl=17   #minimum read length

echo "Creating dir..."
# Create directories if needed
mkdir -p -m 777 /data/work/Projects/BrScRNAseq/fastq_trimmed/$TISSUE
mkdir -p -m 777 /data/work/Projects/BrScRNAseq/fastqc_trimmed/$TISSUE

# get FASTQ files in array
ARRAY=($(ls ${FASTQPATH}/*.fastq.gz | grep R2))
printf '%s\n' "${ARRAY[@]}"


echo "Trimming reads..."
# Trim Adapters
for E in "${ARRAY[@]}"; do
    NAMER2=$(basename "$E" .fastq.gz)
    
    R1FASTQ=$(sed 's/R2/R1/g' <<< $E)
    NAMER1=$(basename "$R1FASTQ" .fastq.gz)
    echo "R1:" $NAMER1
    echo "R2:" $NAMER2

    cutadapt -j $CORES -q $trimQ -A "A{55}" -G file:/data/work/Projects/BrScRNAseq/data/10x5pAdapters.fasta --minimum-length $minl -o /data/work/Projects/BrScRNAseq/fastq_trimmed/${TISSUE}/${NAMER1}.fastq.gz -p /data/work/Projects/BrScRNAseq/fastq_trimmed/${TISSUE}/${NAMER2}.fastq.gz $R1FASTQ $E  > /data/work/Projects/BrScRNAseq/fastqc_trimmed/${TISSUE}/${NAMER2}_cutadapt.log

done

echo "DONE! >:v"



