#!/bin/bash
#SBATCH --partition=genetics_1       # Partition (job queue)
#SBATCH --requeue                    # Return job to the queue if preempted
#SBATCH --job-name=hisat2          # Assign a short name to your job
#SBATCH --nodes=1                    # Number of nodes you require
#SBATCH --ntasks=1                   # Total # of tasks across all nodes
#SBATCH --cpus-per-task=1            # Cores per task (>1 if multithread tasks)
#SBATCH --mem=120G                   # Real memory (RAM) required (MB)
#SBATCH --time=72:00:00              # Total run time limit (HH:MM:SS)
#SBATCH --output=hisat_index.%N.%j.out     # STDOUT output file
#SBATCH --error=hisat_index.%N.%j.err      # STDERR output file (optional)
#SBATCH --export=ALL                 # Export you current env to the job env

#HISAT2...
module load HISAT2/2.1.0

hisat2-build -p 16 DtakHiC1_allchrom.fasta DtakHiC1_allchrom
hisat2-build -p 16 DkikHiC1_allchrom.fasta DkikHiC1_allchrom
hisat2-build -p 16 DanaHiC1_allchrom.fasta DanaHiC1_allchrom
hisat2-build -p 16 DbipHiC1_allchrom.fasta DbipHiC1_allchrom
