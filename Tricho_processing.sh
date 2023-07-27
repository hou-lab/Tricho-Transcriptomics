#!/bin/bash


INTERLEAVE=/home/shengwei/ESP_Metatrans/scripts/merge-paired-reads.sh
UNMERGE=/home/shengwei/ESP_Metatrans/scripts/unmerge-paired-reads.sh
rRNA_databases=/home/shengwei/software/src/sortmerna/sortmerna-2.1b/rRNA_databases
sortmerna=/home/shengwei/software/src/sortmerna/sortmerna-2.1b/sortmerna
diamond_nr_dir=/home/db/blastdb/diamond
diamond_nr=/home/db/blastdb/diamond/nr

uniprot_folder=/home/db/uniprot
Pfam_folder=/home/db/Pfam

# --------
#  FastQC
# --------
DEMUL_FOLDER="raw_reads"
RAW_READS=$DEMUL_FOLDER
FASTQC_FOLDER="02_demultiplexed_FastQC"
mkdir -p $FASTQC_FOLDER

:<<'COMMENT'
source activate NGSprep
for f in `ls $DEMUL_FOLDER/*.fastq.gz`; do
    BASENAME=$(basename $f)
    FILESTEM="${BASENAME%.fastq.gz}"
    if [[ ! -f ${FASTQC_FOLDER}/${FILESTEM}_fastqc.html ]]; then
        fastqc -f fastq -t 5 $f -o $FASTQC_FOLDER # 2>&1 | tee -a ${FASTQC_FOLDER}/${FILESTEM}_FastQC_stdout_stderr.log
    fi
done
source deactivate
# count read number
cd $FASTQC_FOLDER
for f in `ls *.zip`; do
    /home/shengwei/scripts/count_total_reads_number.py $f -t fastqc -c zip -f
done
cd ../
COMMENT


# ----------------------------------
# Adapter and Quality Trimming
# ----------------------------------
TRIMMING_RESULTS="03_adapter_quality_trimming"
mkdir -p $TRIMMING_RESULTS

FASTQC_FOLDER="04_fastqc_after_trimming"
mkdir -p $FASTQC_FOLDER

:<< 'COMMENT'
source activate NGSprep
# only specify the front part of index adapter will work
# 'A' need to be added to the index adapter since it's introduced by complementary 'T' in the Universal adapter
for f in `ls $DEMUL_FOLDER/*.fastq.gz`; do
    BASENAME=$(basename $f)
    FILESTEM="${BASENAME%.fastq.gz}"
    echo "--------"

    atropos \
        -b GGTACACGACGCTCTTCCGATCTGGTCCGG \
        -b TACACGACGCTCTTCCGATCTGGCGCTCATT \
        -b TACACGACGCTCTTCCGATCTGCGAGATTC \
        -b AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -B AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -B AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA \
        -e 0.1 \
        --quality-base 33 \
        --quality-cutoff 20 \
        --trim-n \
        --minimum-length 30 \
        -o $TRIMMING_RESULTS/${FILESTEM}_atropos.fastq \
        --format fastq  --quality-base 33 \
        -se $f \
        --report-file $TRIMMING_RESULTS/${FILESTEM}_atropos_report.txt --report-formats txt \
        --threads 10
    sickle se -f $TRIMMING_RESULTS/${FILESTEM}_atropos.fastq -t sanger \
        -o $TRIMMING_RESULTS/${FILESTEM}_atropos_sickle.fastq \
        --qual-threshold 10 \
        --length-threshold 30 2>&1 | tee -a $TRIMMING_RESULTS/${FILESTEM}_atropos_sickle_stdout_stderr.log

    # run fastqc and gzip fastq
    fastqc -f fastq -t 4 $TRIMMING_RESULTS/${FILESTEM}_atropos_sickle.fastq -o $FASTQC_FOLDER  && gzip $TRIMMING_RESULTS/${FILESTEM}_atropos_sickle.fastq & 
    rm $TRIMMING_RESULTS/${FILESTEM}_atropos.fastq
done

# count read number
cd $FASTQC_FOLDER
for f in `ls *.zip`; do
    /home/shengwei/scripts/count_total_reads_number.py $f -t fastqc -c zip -f
done
cd ../
source deactivate
COMMENT


# ---------------
# Decontamination
# ---------------
# prepare bbmap human ref-genome: http://seqanswers.com/forums/showthread.php?t=42552
# 

BBMAP_RESOURCES="/home/shengwei/software/src/bbmap/resources/"
DECON_FOLDER="05_bbmap_decontamination"
mkdir -p $DECON_FOLDER

:<<'COMMENT'
source activate NGSprep
for f in `ls $TRIMMING_RESULTS/*.fastq.gz`; do
        read=$f
        BASENAME=$(basename $read)
        FILESTEM="${BASENAME%.fastq.gz}"
        echo "--------"
        echo $read

        # remove phiX174 and sequencing artifacts 
        bbduk.sh -Xmx40g threads=10 k=31 \
                 in=$read \
                 ref=${BBMAP_RESOURCES}/phix174_ill.ref.fa.gz,${BBMAP_RESOURCES}/sequencing_artifacts.fa.gz \
                 out=${DECON_FOLDER}/${FILESTEM}_nocon.fq \
                 outm=${DECON_FOLDER}/${FILESTEM}_con.fq \
                 stats=${DECON_FOLDER}/${FILESTEM}_con_stats.txt refstats=${DECON_FOLDER}/${FILESTEM}_con_refstats.txt \
                 statscolumns=5 2>&1 | tee ${DECON_FOLDER}/${FILESTEM}_con_stdout_stderr.log         

        # remove human genome sequences
        bbmap.sh -Xmx40g threads=10 minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast \
                 minhits=2 path=/home/db/bbmap_remove_human_contamination_db/ \
                 qtrim=rl trimq=10 untrim \
                 in=${DECON_FOLDER}/${FILESTEM}_nocon.fq \
                 outu=${DECON_FOLDER}/${FILESTEM}_nocon_nohuman.fq \
                 outm=${DECON_FOLDER}/${FILESTEM}_nocon_human.fq 2>&1 | tee ${DECON_FOLDER}/${FILESTEM}_con_human_stdout_stderr.log

        # clean up
        rm ${DECON_FOLDER}/${FILESTEM}_nocon_human.fq ${DECON_FOLDER}/${FILESTEM}_nocon.fq ${DECON_FOLDER}/${FILESTEM}_con.fq
        gzip ${DECON_FOLDER}/${FILESTEM}_nocon_nohuman.fq & 
done
source deactivate
COMMENT


SSU_PROFILING_RESULTS="06_SSU_profiling_results"
mkdir -p $SSU_PROFILING_RESULTS


# ----------------- #
# run phyloFlash    #
# ----------------- #
phyloFlash_RESULTS=$SSU_PROFILING_RESULTS/"phyloFlash"
mkdir -p ${phyloFlash_RESULTS}

phyloFlash_dbhome="/home/shengwei/Pier_DTS/test_phyloflash/phyloFlash/132/"
phyloFlash_dbhome="/home/db/phyloFlash_Silva_132/"

# change read length parameter
:<<'COMMENT'
source activate phyloFlash
for f in `ls ${DECON_FOLDER}/*_atropos_sickle_nocon_nohuman.fq.gz`; do
    fwd_read=$f
    rev_read=${fwd_read/fwd/rev}
    BASENAME=$(basename $fwd_read)
    FILESTEM="${BASENAME%_atropos_sickle_nocon_nohuman.fq.gz}"
    #FILESTEM="${BASENAME%_nocon_nohuman_fwd.fq}"
    #sampleID=$(echo $FILESTEM | cut -d"_" -f 1,2)
    echo "--------"
    echo $fwd_read
    echo $rev_read

    if [[ $FILESTEM =~ "_non_rRNA_" ]]; then
        echo "non-rRNA reads, continue ..."
        continue
    fi

    read_length=71

    phyloFlash.pl -lib ${FILESTEM} -CPUs 20 -readlength ${read_length} -read1 ${fwd_read} -dbhome ${phyloFlash_dbhome} -treemap -html -skip_emirge 2>&1 | tee ${phyloFlash_RESULTS}/${FILESTEM}_phyloFlash_stdout_stderr.log
    mv ${FILESTEM}.* ${phyloFlash_RESULTS}
    
done
source deactivate
COMMENT

# make heatmap and barplots
:<<'COMMENT'
source activate phyloFlash
cd ${phyloFlash_RESULTS}
phyloFlash_gz_files=`ls -1v *.tar.gz | paste -sd "," -`
for i in $(seq 1 7); do
    phyloFlash_compare.pl --level ${i} --zip ${phyloFlash_gz_files} --task heatmap,barplot,matrix --outfmt pdf --out Tricho_phyloFlash_compare_lvl_${i} --log
done
cd ../../
source deactivate
COMMENT


# ----------
# SortMeRNA
# ----------

sortMeRNA_RESULTS="11_MT_sortMeRNA_results"
mkdir -p $sortMeRNA_RESULTS

:<<'COMMENT'
for f in `ls $DECON_FOLDER/*_atropos_sickle_nocon_nohuman.fq.gz`; do
    BASENAME=$(basename $f)
    FILESTEM="${BASENAME%_atropos_sickle_nocon_nohuman.fq.gz}"
    echo "--------"
    echo $f

    gunzip -ck $f > ${FILESTEM}.fastq

    # sortmerna
    # alternatively, use merged database --ref ${rRNA_databases}/all_rRNAs.fasta,${rRNA_databases}/all_rRNAs.idx \
    $sortmerna --ref ${rRNA_databases}/silva-bac-16s-id90.fasta,${rRNA_databases}/silva-bac-16s-id90.idx:${rRNA_databases}/silva-bac-23s-id98.fasta,${rRNA_databases}/silva-bac-23s-id98.idx:${rRNA_databases}/silva-arc-16s-id95.fasta,${rRNA_databases}/silva-arc-16s-id95.idx:${rRNA_databases}/silva-arc-23s-id98.fasta,${rRNA_databases}/silva-arc-23s-id98.idx:${rRNA_databases}/silva-euk-18s-id95.fasta,${rRNA_databases}/silva-euk-18s-id95.idx:${rRNA_databases}/silva-euk-28s-id98.fasta,${rRNA_databases}/silva-euk-28s-id98.idx:${rRNA_databases}/rfam-5s-database-id98.fasta,${rRNA_databases}/rfam-5s-database-id98.idx:${rRNA_databases}/rfam-5.8s-database-id98.fasta,${rRNA_databases}/rfam-5.8s-database-id98.idx \
               --reads ${FILESTEM}.fastq --num_alignments 1 \
               --fastx --aligned ${sortMeRNA_RESULTS}/${FILESTEM}_rRNA --other ${sortMeRNA_RESULTS}/${FILESTEM}_non_rRNA \
               --log -a 10 -m 24000 -v 2>&1 | tee ${sortMeRNA_RESULTS}/${FILESTEM}_SortMeRNA_stdout_stderr.log
    rm ${FILESTEM}.fastq
done
COMMENT


# ------------ #
# read sam grp #
# ------------ #

SAM_GRP_RESULTS="12_SAM_GRP_results"
mkdir -p $SAM_GRP_RESULTS

reference_fna="./genome/NC_008312.1.fa"
reference_fna="./genome/GCA_000014265.1_ASM1426v1.fna"
reference_gff="./genome/GCA_000014265.1_ASM1426v1.gff"
ID=95

:<<'COMMENT'
source activate py37
for f in `ls ${sortMeRNA_RESULTS}/*_non_rRNA.fastq`; do
    ./generate_sam_grp.py ${f} ${reference_fna} ${ID} -o ${SAM_GRP_RESULTS}
done
conda deactivate
COMMENT


# ------------- #
# featureCounts #
# ------------- #
featureCounts_RESULTS="13_featureCounts_results"
mkdir -p $featureCounts_RESULTS

sam_files=`ls -1v ${SAM_GRP_RESULTS}/after_segemehl/*-rep*_non_rRNA_mapped_to_*.sam | grep -v "PSS"| grep -v "TSS"|  paste -sd " " -`
threads=20

#:<<'COMMENT'
featureCounts -a ${reference_gff} \
                  -o ${featureCounts_RESULTS}/Tricho_featureCounts_result.tab \
                  -t gene \
                  -g locus_tag \
                  -s 1 \
                  -Q 12 \
                  -T $threads \
                  --primary --ignoreDup -d 30 \
		  ${sam_files}
#COMMENT




