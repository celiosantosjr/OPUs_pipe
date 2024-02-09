#!/bin/bash

# download references and expected outputs
OPUs_testingset = 'https://zenodo.org/records/10641930/files/OPUs_testingset.tar.xz?download=1&preview=1'
wget $OPUs_testingset
tar -xvf OPUs_testingset.tar.xz
rm -rf mv OPUs_testingset OPUs_testingset.tar.xz

# run OPUs_pipe
python3 ../main.py -ofolder testresult/ \
                   -pfolder sequences/ \
                   -annofolder sequences/ \
                   -annoxt _FASTA_emapper.annotations.tsv.gz \
                   -minocc 10 \
                   -xt _FASTA_predicted_cds.faa.gz


# Compare two directories
diff -r expected/ testresult/ > output.txt

# Eliminating files that change the order of rows according parallelization
grep -v 'summed_prots.faa.gz\|result_rep_seq.fasta.xz\|OPUs_cluster_relationship.tsv.xz' output.txt > t
mv t output.txt

# Check if there are any differences
if [ -s output.txt ]
then
    # Print the name of the file not matching in both folders
    echo "ERROR - Installation produced unstable version. Files do not match the expected"
    grep -r "Only in folder1" output.txt | cut -d' ' -f4- > mismatch.txt
    cat output.txt
    cat mismatch.txt
    rm -rf mismatch.txt output.txt
else
    # Print success message
    echo "Your installation test was successful"
fi

rm -rf testresult/ sequences/ expected/
