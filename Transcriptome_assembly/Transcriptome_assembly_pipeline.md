# 1. Download RNA reads from SRA database
A text file containing the SRR number of the RNA reads available in the SRA database is required. See "fetchSRA_input.txt" as an example.

`module load sratoolkit/2.8.1`
`org=/core/labs/Wegrzyn/bikash/passiflora/edulis/shortRNA_reads`
`basedir=$org/raw_reads`
`FILENAME="fetchSRA_input.txt"` 
`count=0`
`while read LINE`
`do`
`  	let count++`
`        fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files $LINE`
`        echo "$LINE"`
`done < $basedir/$FILENAME`
`echo -e "\nTotal lines read = $count"`

# 2. Check quality with fastQC

`module load fastqc/0.11.7`
`fastqc -t 6 *.fastq.gz`

# 3. Trim low quality reads with TrimGalore
`module load TrimGalore/0.6.5`
`module load cutadapt/2.7`
`module load fastqc/0.11.7`
`mkdir trimmed`
`trim_galore --phred33 --fastqc --length 50 -o trimmed --cores 10 --paired ERR2040389.1.fastq.gz ERR2040389.2.fastq.gz`
Repeat the command for other libraries too.

# 4. Transcriptome assembly with Trinity
The script is for assembling multiple RNA libraries independenty in an array, so make sure to assign values for array in the job file. E.g: "#SBATCH --array=0-5"

`org=/core/labs/Wegrzyn/bikash/passiflora/edulis/RNA_reads`
`libraries=($(cat $org/trinity/fetchSRA_input.txt)) # for reading SRR numbers`
`trimmeddir=$org/trimmed`
`module load trinity/2.8.5`
`module load bowtie2/2.3.4.3`
`lib=${libraries[$SLURM_ARRAY_TASK_ID]}`
`PE_1=$trimmeddir/trimmed_${lib}_1.fastq`
`PE_2=$trimmeddir/trimmed_${lib}_2.fastq`
`cd $org/trinity`
`out=trinity_${lib}`
`Trinity --seqType fq --left $PE_1 --right $PE_2 --min_contig_length 300 --output $out --full_cleanup --max_memory 100G --CPU 16`

# 5. Renaming header of Trinity output file
For the species with multiple RNA libraries, the header of the contigs in the Trinity output (prefix.Trinity.fasta) needs to be renamed prior combining them into a single file.
This is not required if there's a single RNA library available for the species.

`ls *.fasta > assemblies.txt`
`prefixdir="trinity_prefix"`
`mkdir $prefixdir`

`while read assembly`
`do`
` 	lib=${assembly#trinity_}`
`        lib=${lib%.Trinity.fasta}`
`        sed "s/>/>${lib}_/g" $assembly > $prefixdir/${lib}.Trinity.fasta`
`done < assemblies.txt`
`rm assemblies.txt`

After renaming the header, concatenate contigs generated from multiple libraries to a single file

`cat $prefixdir/*.Trinity.fasta > $prefixdir/edulis.trinity.fasta`

# 6. Frame selection using EvidentialGene

`module load perl/5.28.1`
`module load blast/2.10.1`
`module load exonerate/2.4.0`
`evigene=/isg/shared/apps/evigene/20190101`
`export PATH=$evigene/scripts/:$evigene/scripts/prot/:$evigene/scripts/rnaseq/:/isg/shared/apps/cdhit/4.6.8/:$PATH`
`org=/core/labs/Wegrzyn/bikash/passiflora/RNA-Seq/preprocessing/PED/trinity`
`tr2aacds.pl -tidy -NCPU 16 -MAXMEM 150000 -log -cdna $org/edulis.trinity.fasta`

# 7. Size selection of ORFs
The open reading frames (ORFs) predicted with the EvidentialGene were further processed to remove sequences less than 300 bp.

`module load seqtk`
`dir1=/core/labs/Wegrzyn/bikash/passiflora/RNA-Seq/preprocessing/PED/evigene/okayset`
`seqtk seq -L 300 "$dir1/edulis.trinity.okay.trinity.cds" > ped_300.cds`
`grep -h ">" ped_300.cds | sed 's/^>//g' > headers_ped_cds.txt`
`seqtk subseq "$dir1/edulis.trinity.okay.trinity.aa" headers_ped_cds.txt > ped_300.pep`

# 8. Transcriptome completeness assessment using BUSCO5
`module load busco/5.0.0`
`module load hmmer/3.2.1`
`module load blast`
`busco -c 12 -m protein -i ped_300.pep -o busco_O -l eudicots_odb10`

# 9. Clustering to remove redundant ORFs
Removal of ORFs that are nearly identical (generated during Trinity assembly as "isoforms") can be performed with clustering with 95% sequence similarity.

`module load vsearch/2.4.3`
`module load seqtk`
`vsearch --threads 8 --log 95LOGFile --cluster_fast ../ped_300.cds --id 0.95 --centroids ped.cds --uc 95clusters.uc`

Then, extract the protein files using non-redundant transcripts generated with clustering        
`grep -h ">" ped.cds | sed 's/^>//g' > header.txt`
`seqtk subseq ped_300.pep header.txt > ped.pep`
`rm header.txt`

# 10. Transcriptome completeness re-assessment using BUSCO5
`module load busco/5.0.0`
`module load hmmer/3.2.1`
`module load blast`
`busco -c 12 -m protein -i ped.pep -o busco_O -l eudicots_odb10`