#!/bin/bash

# Set variables:

# VF_REF=/media/data1/janiek/ref/vf_refseq/GCF_000011805.1_ASM1180v1_genomic.fna
MS2_REF=/media/data1/janiek/ref/MS2_refseq/GCF_000847485.1_ViralProj14659_genomic.fna
# QBETA_REF=/media/data1/janiek/ref/qbeta_knudsen.fasta
R17_REF=/media/data1/janiek/ref/R17.fasta
FR_REF=/media/data1/janiek/ref/fr.fasta
SYNTHETIC_STANDARD_REF=/media/data1/janiek/ref/MS2_synthetic_standard.fasta

FASTQ_DIR=/media/data4to7/janiek/260716_barcoded_RNA_phage_cont/demux_out/unknown/20260717_0436_0_PBM05491_8e943acf/fastq_pass
OUTPUT_DIR=/media/data4to7/janiek/260716_barcoded_RNA_phage_cont/demux_out/unknown/20260717_0436_0_PBM05491_8e943acf/mapping
mkdir -p "$OUTPUT_DIR"

map_one () { 
  minimap2 -ax map-ont -t 10 "$2" "$1" > "$3.sam"
  samtools sort -@ 2 "$3.sam" -o "$3.bam"
  samtools flagstat "$3.bam" > "${3}_flagstat.txt"
  samtools index "$3.bam"
}

# assign the right ref for each barcode
get_ref () {
	local n=$1
	if   (( n >= 1  && n <= 6  )); then echo "$MS2_REF"
	elif (( n >= 7  && n <= 8  )); then echo "$FR_REF"
	elif (( n >= 9  && n <= 10 )); then echo "$R17_REF"
	elif (( n >= 11 && n <= 16 )); then echo "$SYNTHETIC_STANDARD_REF"
	else echo ""      # (no reference assigned)
	fi
}

# mapping
## map barcodes as below: 
## 1-6: MS2_REF
## 7-8: FR_REF
## 9-10: R17_REF
## 11-16: SYNTHETIC_STANDARD_REF

# loop over all barcode directories present
	for BC_DIR in "$FASTQ_DIR"/barcode*; do
	BARCODE=$(basename "$BC_DIR")
	NUM=$((10#${BARCODE#barcode}))

	REF=$(get_ref "$NUM")
	if [[ -z "$REF" ]]; then
	echo "Skipping $BARCODE (no reference assigned)"
	continue
	fi

	echo "Processing $BARCODE -> $(basename "$REF")"

	map_one "$BC_DIR"/*.fastq "$REF" "$OUTPUT_DIR/${BARCODE}"

done

echo "ALL DONE!"
