# 1. Check read quality with fastQC
## Running FastQC
`module load fastqc/0.11.7`
`fastqc -t 6 *.fq.gz`

## Trim low quality reads with TrimGalore
`module load TrimGalore/0.6.5`
`module load cutadapt/2.7`
`module load fastqc/0.11.7`
`mkdir trimmed`
`trim_galore --phred33 --fastqc --length 50 -o trimmed --cores 10 --paired POT_R1.fq.gz POT_R2.fq.gz`

*repeat the commands for other species as well

# 2. Genome Assembly
## Running MaSuRCA

`module load singularity/3.1.1`
`module load perl/5.28.1`
`module load MaSuRCA/4.0.3`

`masurca config_file_pot` #make sure to copy the config file in the working directory and adjust parameters accordingly
`./assemble.sh`

# 3. Genome Annotation

## Removal of contigs < 1500 bp from the MaSuRCA assembled genome
`module load seqtk`
`seqtk seq -L 1500 /core/labs/Wegrzyn/bikash/passiflora/NovaSeq/masurca/masurca_pot/CA/primary.genome.scf.fasta.masked > pot_filtered.fasta.masked`


## Model repeats in the assembled genome using RepeatModeler2
`module load perl/5.24.0`
`module load RepeatModeler/2.01`
`module load genometools/1.6.1`
`module load mafft/7.471`
`module load cdhit/4.8.1`
`module load ninja/0.95`

`export PERL5LIB=/UCHC/PublicShare/szaman/perl5/lib/perl5/`
`export LTRRETRIEVER_PATH=/core/labs/Wegrzyn/annotationtool/software/LTR_retriever`

`org=/core/labs/Wegrzyn/bikash/passiflora/NovaSeq/masurca/masurca_pot/CA`

###BuildDatabase -name POT $org/primary.genome.scf.fasta
`RepeatModeler -database POT -pa 20 -LTRStruct`

## Mask identified repeats using RepeatMasker4

`export PATH=/core/labs/Wegrzyn/annotationtool/software/RepeatMasker/4.0.6:$PATH`
`org=/core/labs/Wegrzyn/bikash/passiflora/NovaSeq/masurca/masurca_pot/CA`

`RepeatMasker -lib POT-families.fa -pa 32 -gff -a -noisy -low -xsmall $org/primary.genome.scf.fasta`

## Annotation of protein coding genes using BRAKER2 (For three species- P. auriculata, P. oerstedii and P. pittieri)
`module load python/3.6.3`
`module load biopython/1.70`
`module load bamtools/2.5.1`
`module load blast/2.10.0`
`module load genomethreader/1.7.1`
`module load perl/5.28.1`

### This script expects these programs to be in these specified folders
### If they are ever moved/updated, or if you want to use ones in your home directory update the paths below
`augustusPath="/core/labs/Wegrzyn/annotationtool/software/Augustus_3.4.0"`
`brakerPath="/core/labs/Wegrzyn/annotationtool/software/BRAKER_2.1.5/scripts"`
`cbdfastaPath="/core/labs/Wegrzyn/annotationtool/software/cdbfasta"`
`genemarkPath="/core/labs/Wegrzyn/annotationtool/software/gmes_linux_64"`

`export PATH=$brakerPath:$brakerPath/scripts:$PATH`
`export CDBTOOLS_PATH=$cbdfastaPath`
`export GENEMARK_PATH=$genemarkPath`
`export BAMTOOLS_PATH=/isg/shared/apps/bamtools/2.5.1/bin`
`export BLAST_PATH=/isg/shared/apps/blast/ncbi-blast-2.10.0+/bin`
`export SAMTOOLS_PATH=/isg/shared/apps/samtools/1.9/bin`
`export BLAST_PATH=/isg/shared/apps/blast/ncbi-blast-2.10.0+/bin`
`export ALIGNMENT_TOOL_PATH=/isg/shared/apps/gth/1.7.3/bin`

`export AUGUSTUS_BIN_PATH=$augustusPath/bin`
`export AUGUSTUS_SCRIPTS_PATH=$augustusPath/scripts`
`export AUGUSTUS_CONFIG_PATH=$HOME/augustus/config`

`speciesname="pot"`
`genome="/core/labs/Wegrzyn/bikash/passiflora/NovaSeq/annotation/masked_filtered/pot_filtered.fasta.masked"`
`bamfile="/core/labs/Wegrzyn/bikash/passiflora/NovaSeq/annotation/POT/hisat2/sorted_pot.bam"`

`mkdir -p tmp`
`export TMPDIR=/core/labs/Wegrzyn/bikash/passiflora/NovaSeq/annotation/PPT/braker/tmp`

`script="$brakerPath/braker.pl --species=$speciesname --genome=$genome --bam=$bamfile --cores 20 --softmasking 1 --gff3"`
`eval $script`

## Generate RNA-alignment to feed as an external evidence for the three species

### Index the soft-masked genome first (For three species- P. auriculata, P. oerstedii and P. pittieri)
`module load hisat2/2.1.0`
`hisat2-build -p 16 /core/labs/Wegrzyn/bikash/passiflora/NovaSeq/annotation/masked_filtered/pot_filtered.fasta.masked pot_masked`

### Align RNA reads against the masked genome
`module load hisat2/2.1.0`
`module load samtools`

`dir=/core/labs/Wegrzyn/bikash/passiflora/RNA-Seq`
`hisat2 -x pot_masked -1 "$dir/pot_R1.fq.gz" -2 "$dir/pot_R2.fq.gz" -p 16 -S pot.sam`
`samtools view -@ 16 -uhS pot.sam | samtools sort -@ 16 -T pot_temp -o sorted_pot.bam`

## Genome completeness assessment using BUSCO5

`module load busco/5.0.0`
`module load hmmer/3.2.1`
`module load blast`

`infile=/core/labs/Wegrzyn/bikash/passiflora/NovaSeq/annotation/POT/braker/braker/augustus.hints.aa`
`outfold=pot_agustus_hints`
`busco -c 16 -m prot -i "$infile" -o "$outfold" -l eudicots_odb10`

