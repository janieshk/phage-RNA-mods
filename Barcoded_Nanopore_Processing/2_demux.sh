#!/bin/bash

/media/data0/janiek/programs/dorado-2.0.1-linux-x64/bin/dorado demux \
	--kit-name SQK-DRB004-24 \
	--emit-fastq \
	--output-dir /media/data4to7/janiek/260716_barcoded_RNA_phage_cont/demux_out/ \
	--threads 20 \
	/media/data4to7/janiek/260716_barcoded_RNA_phage_cont/260716_barcoded_RNA_phage_cont_rebasecall_260618_hac_v3.0.1.bam
