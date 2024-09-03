#!/bin/bash


# ------------
# circolize DE 
# ------------
#./circolize_DE.py ./Tricho_HFe_vs_LFe_contrast_lvl_anno.tsv -n log2FoldChange padj baseMean_HFe baseMean_LFe --geneID_col_name locus_tag
#cat Tricho_HFe_vs_LFe_contrast_lvl_anno_DE_plus_up.txt Tricho_HFe_vs_LFe_contrast_lvl_anno_DE_minus_up.txt > Tricho_HFe_vs_LFe_contrast_lvl_anno_DE_up.txt
#cat Tricho_HFe_vs_LFe_contrast_lvl_anno_DE_plus_down.txt Tricho_HFe_vs_LFe_contrast_lvl_anno_DE_minus_down.txt > Tricho_HFe_vs_LFe_contrast_lvl_anno_DE_down.txt
#cat Tricho_HFe_vs_LFe_contrast_lvl_anno_plus_log2FoldChange.txt Tricho_HFe_vs_LFe_contrast_lvl_anno_minus_log2FoldChange.txt > Tricho_HFe_vs_LFe_contrast_lvl_anno_log2FoldChange.txt
cat Tricho_HFe_vs_LFe_contrast_lvl_anno_plus_baseMean_HFe.txt Tricho_HFe_vs_LFe_contrast_lvl_anno_minus_baseMean_HFe.txt > Tricho_HFe_vs_LFe_contrast_lvl_anno_baseMean_HFe.txt
cat Tricho_HFe_vs_LFe_contrast_lvl_anno_plus_baseMean_LFe.txt Tricho_HFe_vs_LFe_contrast_lvl_anno_minus_baseMean_LFe.txt > Tricho_HFe_vs_LFe_contrast_lvl_anno_baseMean_LFe.txt

# -------------
# circolize fna
# -------------
#./circolize_fna.py GCA_000014265.1_ASM1426v1_genomic.fna -f  

# -------------
# circolize gff
# -------------
#./circolize_gff.py GCA_000014265.1_ASM1426v1_genomic.gff -f


# -------------
# circolize grp
# -------------

:<<'COMMENT'
source activate py37
python ./circolize_grp.py -w 1000 -p HFe-d5_norm_merge \
    HFe-d5_norm_merge.grp CP000393.1
python ./circolize_grp.py -w 1000 -p LFe-d5d7_norm_merge \
    LFe-d5d7_norm_merge.grp CP000393.1
conda deactivate
COMMENT



