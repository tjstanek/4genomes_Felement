#!/bin/bash

#SBATCH --partition=main       # Partition (job queue)
#SBATCH --requeue                     # Return job to the queue if preempted
#SBATCH --job-name=orthofind             # Assign an short name to your job
#SBATCH --nodes=1                     # Number of nodes you require
#SBATCH --ntasks=1                    # Total # of tasks across all nodes
#SBATCH --cpus-per-task=12            # Cores per task (>1 if multithread tasks)
#SBATCH --mem=120G                    # Real memory (RAM) required (MB)

#SBATCH --time=72:00:00               # Total run time limit (HH:MM:SS)
#SBATCH --output=orthofinder.%N.%j.out  # STDOUT output file
#SBATCH --error=orthofinder.%N.%j.err   # STDERR output file (optional)
#SBATCH --export=ALL                  # Export you current env to the job env

orthofinder -f 4genomes_faa
