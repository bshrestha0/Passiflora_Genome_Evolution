# 1. Identification of single-copy orthologues using OrthoFinder
The input for the OrthoFinder run is a directory containing amino acid sequences for the non-redundant protein-coding genes for all the species.

`module load OrthoFinder/2.5.1`
`module load DLCpar/1.0`
`module load FastME`
`module load diamond/0.9.25`
`module load mcl`
`module load mafft/7.471`
`module load fasttree/2.1.10`
`org=/core/labs/Wegrzyn/EASEL/passiflora/analysis/orthofinder`
`orthofinder  -f $org/proteinFiles -S diamond -t 12`

The single copy orthologues will in the the directory "Single_Copy_Orthologue_Sequences" within the OrthoFinder output directory.

# 2. File preparation for substitution rate analyses
## Generating protein and CDS alignments for the orthogroups.

First align the single copy orthologue sequences for each orthogroup using MAFFT.

`module load mafft/7.471
`org=/core/labs/Wegrzyn/EASEL/passiflora/analysis/orthofinder/OrthoFinder_output/Single_Copy_Orthologue_Sequences`
`outDir=/core/labs/Wegrzyn/EASEL/passiflora/analysis/single_copy_rate/alignment`
`for file in $org/*fa
do
name=${file%.*}
mafft --auto $file > $outDir/"$name".pep.fasta
done`
It will generate alignment for all the orthogroups with filenames ending with ".pep.fasta".

Collecting CDS from the corresponding protein orthogroups alignments. All the CDS for the protein sequences can be concatenated to generate a single file.
`cat *.cds > CDS`

Extract the sequence header (ID) from the protein alignment files using for loop.
`
for file in *.fasta
do
name=${file%.*}
grep -e ">" $file > $name.txt
done
`
Use the file with header to extract CDS per orthogroup from the single concatenated CDS file.
`
for header in *.txt
do
part=${header%.*}
while read p;
do  grep -A1 "^$p$" CDS >> $part.cds ;
done < "$header"
done
`
Use PAL2NAL to convert multiple sequence alignment of proteins and DNA sequence to generate codon alignment.
PAL2NAL can be downloaded from http://www.bork.embl.de/pal2nal/index.cgi?example=Yes#Download
`
for file in *.pep.fasta
do
part=${file%%.*}
/home/FCAM/bshrestha/pal2nal.v14/pal2nal.pl  $part.pep.fasta $part.cds -output paml > $part.cds.fasta
done
`

## Filtering orthogroups for substitution rate analysis

Remove the orthogroups if it contains sequence(s) that are shorter than 70% of the total alignment length. These orthogroups will be discarded from the PAML analysis.
Run "fasta_drop.py" to generate new alignments that lacks the sequences shorter thatn 70% of the total alignment length.

`module load python/3.8.1`
`module load biopython/1.70`
`org=/core/labs/Wegrzyn/EASEL/passiflora/analysis/single_copy_rate/script`
`
for file in *fasta
do
name=${file%%.*}
python3.8 $org/fasta_drop.py $file $name.drop 0.7
done
`
If the number of sequences in the new alignment is less than the original alignment, remove those orthogroups.

`
for file in *fasta
do
name=${file%%.*}
ori=$(grep -c ">" $file)
edit=$(grep -c ">" $name.drop)
  if [ "$ori" != "$edit" ]
  then
    echo $file
    rm $file
fi
done
`

# 3. Pairwise substitution rate analysis using PAML

Pairwise substitution rate analysis requires alignment file for each orthogroup and the parameter file for PAML run
First, generate a file (ID.txt) that contains orthogroup name, which can be used for generating individual directory for the pairwise run.
 
`for file in *.cds.fasta
do
OG=${file%%.*}
cat $OG > ID.txt
done`

Use the orthogroup ID to create a directory with the alignment and mlc (parameter) file for PAML.

`while read line
do
mkdir $line
cp codeml.ctl "$line".ctl
sed -i "s/input_fasta/${line}.cds.fasta/g" "$line".ctl #rename fasta file
sed -i "s/mlc_file/mlc_${line}/g" "$line".ctl #rename mlc filename
mv ${line}.ctl ${line}
cp /core/labs/Wegrzyn/EASEL/passiflora/analysis/single_copy_rate/final/${line}.* ${line}
done < ID.txt 
`
Running PAML for the pairwise analysis

`module load  paml/4.9`
`while read id
do
cd ${id}
codeml ${id}\.ctl
cd ../
done < ID.txt`

# 3. Branch-specific (free-ratio) substitution rate analysis using PAML

Branch-specific substitution rate analysis requires a tree file in addition to alignment and parameter file.
I am using Maximum Likelihood tree (5genes.tre) generated with IQ-TREE using five low-copy nuclear genes.
File preparation for the free-ratio run includes creating directory with alignment file, parameter file and the tree file.

`
while read line
do
mkdir $line
cp codeml.ctl "$line".ctl
sed -i "s/input_fasta/${line}.cds.fasta/g" "$line".ctl #rename fasta file
sed -i "s/mlc_file/mlc_${line}/g" "$line".ctl #rename mlc filename
#copying files
mv ${line}.ctl ${line}
cp /core/labs/Wegrzyn/EASEL/passiflora/analysis/single_copy_rate/final/${line}.* ${line}/
cp 5genes.tre ${line}/
#rename the sequences in the fasta so that it match sequence name in the tree file.
sed -i "s/aha.*/aha/g" ${line}/${line}.cds.fasta
sed -i "s/AT[0-9]G.*/ath/g" ${line}/${line}.cds.fasta
sed -i "s/dco.*/dco/g" ${line}/${line}.cds.fasta
sed -i "s/dter.*/dte/g" ${line}/${line}.cds.fasta
sed -i "s/Mane.*/mea/g" ${line}/${line}.cds.fasta
sed -i "s/pas.*/pas/g" ${line}/${line}.cds.fasta
sed -i "s/pba.*/pba/g" ${line}/${line}.cds.fasta
sed -i "s/ped.*/ped/g" ${line}/${line}.cds.fasta
sed -i "s/pmo.*/pmo/g" ${line}/${line}.cds.fasta
sed -i "s/pms.*/pms/g" ${line}/${line}.cds.fasta
sed -i "s/pos.*/pos/g" ${line}/${line}.cds.fasta
sed -i "s/Potr.*/pta/g" ${line}/${line}.cds.fasta
sed -i "s/ppt.*/ppt/g" ${line}/${line}.cds.fasta
sed -i "s/pot\..*/pot/g" ${line}/${line}.cds.fasta
sed -i "s/tjo.*/tjo/g" ${line}/${line}.cds.fasta
sed -i "s/tsu.*/tsu/g" ${line}/${line}.cds.fasta
done < ID.txt

Running PAML for the free-ratio analysis

`module load  paml/4.9`
`while read id
do
cd ${id}
codeml ${id}\.ctl
cd ../
done < ID.txt

# 4. Functional annotation of single copy orthologues

Using the ID of the single copy orthologues, extract protein sequences for Arabidopsis thaliana from each orthogroup.

`org=/core/labs/Wegrzyn/EASEL/passiflora/analysis/orthofinder/OrthoFinder_output/Single_Copy_Orthologue_Sequences`
`
while read line
do
grep -A0 ">AT" ${org}/${line}\.fa > header
sed "s/$/_${line}/g" header >> ATHsingleCopy.fa
grep -A1 ">AT" ${org}/${line}\.fa | grep -v ">AT" >> ATHsingleCopy.fa
done < ID.txt
`

Annotate using Entap

`module load anaconda/2.4.0`
`module load perl/5.24.0`
`module load diamond/0.9.36`
`module load interproscan/5.25-64.0`
`prot=/core/labs/Wegrzyn/EASEL/passiflora/analysis/single_copy_rate/single_copy/ATHsingleCopy.fa`

`/core/labs/Wegrzyn/EnTAP/EnTAP_v0.10.4/EnTAP/EnTAP --runP \
--ini entap_config.ini \
-i $prot \
-d /isg/shared/databases/Diamond/Uniprot/uniprot_sprot.dmnd \
-d /isg/shared/databases/Diamond/RefSeq/plant.protein.faa.205.dmnd \
--threads 16 
`
`
