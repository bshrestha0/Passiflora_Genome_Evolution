# Estimation of whole genome duplication based on synonymous substitution distribution of protein coding genes (https://github.com/arzwa/wgd).
For the analysis, an entire set of non-redundant protein coding genes is required for each species.

WGD software (Zwaenepoel and Van de Peer, 2019) was installed using Anaconda v.3

`module load anaconda3/4.1.1`
`source activate wgd_env`
`module unload anaconda3/4.1.1`

`module load mcl/14-137 `
`module load blast/2.11.0`
`module load fasttree/2.1.10`
`module load mafft/7.471`
`module load paml/4.9`

`data=/core/labs/Wegrzyn/bikash/passiflora/analysis/gene_set`

All-vs-all blast using blastP and MCL clustering. The E-value cut-off was set to 1e−10 with an inflation factor of 2.0.
`wgd mcl --cds --mcl -s $data/aha.cds  -o ./wgd_blast -n 8`

Ks distribution using fastTree node weighting
`mkdir Syn`
`wgd ksd wgd_blast/aha.cds.blast.tsv.mcl $data/aha.cds -o ./Syn -n 8`

Fitting Gaussian mixed model (GMM) and Bayesian Gaussian mixed model BGMM
`wgd mix Syn/aha.cds.ks.tsv -n 1 5`
`wgd mix --method bgmm Syn/aha.cds.ks.tsv`

