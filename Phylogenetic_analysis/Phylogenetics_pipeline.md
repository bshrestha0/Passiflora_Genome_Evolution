# 1. Steps for mining low-copy nuclear genes (MSH1, MCM5, MLH1, SMC1 and SMC2)
First download the protein sequences for Arabidopsis thaliana using https://www.arabidopsis.org/tools/bulk/protein/index.jsp 
following the accession numbers: AT2G07690.1, AT3G24320.1, AT4G09140.1, AT3G54670.1, and AT5G62410.1.
Then concatenate all the sequences into a single file, let's say "5genes.pep"". Run the following commands to mine the orthologs in other species.

`module load blast/2.2.29`
`module load samtools/1.9`
`db=/core/labs/Wegrzyn/bikash/passiflora/analysis/phylogeny/dataset/5genes`

`makeblastdb  -in 59genes.pep -out ath_pep -dbtype prot` #make blast database for 5 protein sequences
`blastp -query ./PED/ped.evigene.pep -db $db/ath_pep -evalue 1e-5 -num_threads 6 -max_target_seqs 1 -outfmt 6 -max_hsps 1 -out ped.blastp`

`awk  '$3 > 50 {print $0}' ped.blastp > ped.hits #retain ones with percentage id > 50`
`awk  '{print $1}' ped.hits | xargs samtools faidx ped.cds > ped.5genes.fasta` #use the ID to extract DNA sequences

# 2. Align each gene individually using MAFFT
`module load mafft/7.471`
`mafft --auto MSH1.cds > MSHA.align.cds` #MSH1.cds represents the coding sequences for MSH1 gene for all the species
Repeat the step for other genes as well and concatenate all the alignment into single multiple sequence alignment file.

# 3. Phylogenetic inference

`module load iqtree/2.1.3`
`iqtree2 -s 5genes.phy -p 5genes_partition.nex -st DNA -m MFP -b 200 -nt AUTO -seed 142241 --sampling GENESITE -ninit 2000 -ntop 400 -nbest 100 -nstop 200`
The data is partition by gene and codon positions (5genes_partition.nex).
