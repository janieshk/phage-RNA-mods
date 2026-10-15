#!/bin/bash
#SBATCH --job-name=rebasecall_260718
#SBATCH --output=job_rebasecall_260718.out
#SBATCH --error=job_rebasecall_260718.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --gpus=1
#SBATCH --partition=gpu
#SBATCH -C GPU_SKU:H100_SXM5

module load python/3.9

/home/groups/nicolemm/janiek/dorado-0.9.1-linux-x64/bin/dorado basecaller /home/groups/nicolemm/janiek/dorado-0.9.1-linux-x64/bin/models/rna004_130bps_hac@v3.0.1 \
/oak/stanford/groups/nicolemm/janiek/260716_barcoded_RNA_phages_cont/no_sample_id/20260716_2134_P2S-01916-B_PBM05491_208dca83/pod5 \
--emit-moves \
> /oak/stanford/groups/nicolemm/janiek/260716_barcoded_RNA_phages_cont/no_sample_id/20260716_2134_P2S-01916-B_PBM05491_208dca83/analysis/260716_barcoded_RNA_phage_cont_rebasecall_260618_hac_v3.0.1.bam
