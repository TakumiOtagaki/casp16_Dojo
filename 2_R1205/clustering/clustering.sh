. ~/.bashrc
ml rosetta

cd /large/otgk/casp/casp16/2_R1205/

rm clustering/RNA_clustered.out
rna_cluster.default.linuxgccrelease \
    -in:file:silent R1205.out \
    -cluster:radius 2.0 \
    -out:file:silent clustering/RNA_clustered.out \
    -nstruct 10 # クラスタ数

echo "silent file to pdb"
mkdir cluster_models
cd cluster_models
rna_extract.linuxgccrelease \
 -in:file:silent ../clustering/RNA_clustered.out \
 -in:file:silent_struct_type rna \