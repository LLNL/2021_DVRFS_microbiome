#  03. Phylogenetic Tree

wd=<path to working dir>
cd ${wd}

## Make fasta file with sequences of OTUs that are being kept
input=${wd}/tax_table_for_tree.txt #This is an excel file with OTU ID, sequences, and taxonomy
output=${wd}/sequences_for_tree_onlyOTU.fa
python3 tab2fasta.py ${input} 2 1 > ${output}
sed -e '1,2d' ${output} > ${output}2 #delete first two lines from file
rm ${output}
mv ${output}2 ${output}

## SINA alignment (v1.6.0)
OTU_input=${wd}/sequences_for_tree_onlyOTU.fa
concate_align=${OTU_input%.*}.clean.fa
ref_db=<path to SILVA_138_SSURef_NR99_05_01_20_opt.arb>

conda activate sina

sina -i ${concate_align} -o ${concate_align%.*}_SINAaligned.fasta -r ${ref_db} \
-v --log-file ${concate_align%.*}_SINAaligned_log.txt --fasta-write-dna \
--search --meta-fmt=csv \
--search-max-result 1 --lca-fields=tax_slv

## remove gaps in columns containing only gaps
trimal -in ${concate_align%.*}_SINAaligned.fasta -out ${concate_align%.*}_SINAaligned_removedgaps.fasta -htmlout ${concate_align%.*}_SINAaligned_removedgaps.html -fasta -noallgaps

sed '/^>/ s/ .*//' ${concate_align%.*}_SINAaligned_removedgaps.fasta > ${concate_align%.*}_SINAaligned_removedgaps2.fasta
mv ${concate_align%.*}_SINAaligned_removedgaps2.fasta ${concate_align%.*}_SINAaligned_removedgaps.fasta

## FastTree for checking
fasttree -nt -gtr -log ${concate_align%.*}_SINAaligned_removedgaps_tree.log ${concate_align%.*}_SINAaligned_removedgaps.fasta > ${concate_align%.*}_SINAaligned_removedgaps_tree.tre

## Upload to CIPRES Science Gateway for RAxML Tree
