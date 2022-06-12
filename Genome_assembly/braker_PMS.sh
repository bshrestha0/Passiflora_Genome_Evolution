#!/bin/bash
#SBATCH --job-name=braker_pms
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 20
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=200G
#SBATCH --mail-user=bikash.shrestha@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

module load python/3.6.3
module load biopython/1.70
module load bamtools/2.5.1
module load blast/2.10.0
module load genomethreader/1.7.1
module unload perl
module load perl/5.28.1

## This script expects these programs to be in these specified folders
## If they are ever moved/updated, or if you want to use ones in your home directory update the paths below
augustusPath="/core/labs/Wegrzyn/annotationtool/software/Augustus_3.4.0"
brakerPath="~/BRAKER-2.1.6/scripts"
cbdfastaPath="/core/labs/Wegrzyn/annotationtool/software/cdbfasta"
genemarkPath="/core/labs/Wegrzyn/annotationtool/software/gmes_linux_64"
prothintpath="/home/FCAM/abhattarai/ProtHint/bin"

export PATH=$brakerPath:$brakerPath/scripts:$PATH
export DIAMOND_PATH=/home/FCAM/abhattarai/0.9.25/bin
export CDBTOOLS_PATH=$cbdfastaPath
export GENEMARK_PATH=$genemarkPath
export BAMTOOLS_PATH=/isg/shared/apps/bamtools/2.5.1/bin
export BLAST_PATH=/isg/shared/apps/blast/ncbi-blast-2.10.0+/bin
export SAMTOOLS_PATH=/isg/shared/apps/samtools/1.9/bin
export BLAST_PATH=/isg/shared/apps/blast/ncbi-blast-2.10.0+/bin
export ALIGNMENT_TOOL_PATH=$prothintpath

export AUGUSTUS_BIN_PATH=$augustusPath/bin
export AUGUSTUS_SCRIPTS_PATH=$augustusPath/scripts
export AUGUSTUS_CONFIG_PATH=$HOME/augustus/config

export PROTHINT_PATH=$prothintpath
speciesname="PMS_TEST"
genome="/core/labs/Wegrzyn/bikash/passiflora/NovaSeq/annotation/masked_filtered/pms_filtered.fasta.masked"
protfile="/core/labs/Wegrzyn/bikash/passiflora/NovaSeq/annotation/PMS/getting_prots_for_braker/odb_prot_with_pass_prot.fasta"
prothints="/core/labs/Wegrzyn/bikash/passiflora/NovaSeq/annotation/PMS/braker/braker_worked/prothint_augustus.gff"

mkdir -p tmp
export TMPDIR=/core/labs/Wegrzyn/bikash/passiflora/NovaSeq/annotation/PMS/braker/tmp

script="$brakerPath/braker.pl --species=$speciesname --genome=$genome --hints=$prothints --cores 20 --softmasking --gff3"
eval $script