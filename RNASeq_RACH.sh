#!/bin/bash/

#Usage: sh RNASeq_RACH.sh path/to/STARgenomeindexfolder path/to/gtfAnnotationfile

#Create environment & install Trimmomatic (Anaconda), STAR (https://github.com/alexdobin/STAR) & StringTie (Anaconda)
#Make sure genome is indexed using STAR 
#Only include samples which passed QC run 

cores=18
genome_index=$1
annotation=$2

#Define directory paths, change accordingly
INPUT_DIR="/full/path/to/input/folder"
OUTPUT_DIR="/full/path/to/output/folder"
ADAPTERS="/full/path/to/folder/for/adapterfile"
TRIMMOMATIC="/full/path/to/bin/trimmomatic"
STAR="/full/path/to/STAR-2.7.11a/source/STAR"
STRINGTIE="/full/path/to/bin/stringtie"

#Create output folders (-p, only if not already made)
mkdir -p trimmed_reads/
mkdir -p star_alignment/
mkdir -p stringtie_expression/

#Create function to check success of previous command

check_success() {
  if [ $? -ne 0 ]; then
    echo "Error: The previous command failed. Exiting..."
    exit 1
  fi
}

#Step 1: Quality trim reads using Trimmomatic v0.39

for file in $INPUT_DIR/*_1.fq.gz; do 
     base=$(basename $file _1.fq.gz)
     echo "starting trimming for $base"
     $TRIMMOMATIC PE -threads $cores $INPUT_DIR/${base}_1.fq.gz $INPUT_DIR/${base}_2.fq.gz \
	-baseout $OUTPUT_DIR/trimmed_reads/${base}_trimmed.fq.gz \
	ILLUMINACLIP:$ADAPTERS/all_truseq_edited.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
     check_success

done

#Step 2 :Alignment using STAR v2.7.11a
for file in $OUTPUT_DIR/trimmed_reads/*trimmed_1P.fq.gz; do
    base=$(basename ${file} _trimmed_1P.fq.gz)
    echo "starting alignment for $base"
    $STAR --runThreadN $cores \
        --genomeDir $genome_index \
        --readFilesIn $OUTPUT_DIR/trimmed_reads/${base}_trimmed_1P.fq.gz $OUTPUT_DIR/trimmed_reads/${base}_trimmed_2P.fq.gz \
        --outFileNamePrefix $OUTPUT_DIR/star_alignment/${base}_ \
	--readFilesCommand zcat \
	--quantMode GeneCounts \
        --outSAMtype BAM SortedByCoordinate
    check_success
done

#Step 3:Transcript assembly with StringTie v2.2.1 
for file in $OUTPUT_DIR/star_alignment/*_Aligned.sortedByCoord.out.bam; do
    base=$(basename ${file} _Aligned.sortedByCoord.out.bam)
    echo "starting transcript assembly for $base"
    $STRINGTIE $file -p $cores -G $annotation -e -o $OUTPUT_DIR/stringtie_expression/${base}.gtf\
     -A $OUTPUT_DIR/stringtie_expression/${base}.gene.abundances.tsv -B
    check_success
done

echo "RNA-Seq analysis pipeline completed!"
