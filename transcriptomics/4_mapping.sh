#!/bin/bash
exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>mapping.log 2>&1

CORENUM=60
MAIN_DIR=$(python -c "print('/'.join('${BASH_SOURCE[0]}'.split('/')[:-2]))")
FILES_DIR=$MAIN_DIR/files
INDEX_DIR=$MAIN_DIR/indexes
MAPPING_DIR=$MAIN_DIR/mappings
READS_DIR=$MAIN_DIR/reads

mkdir $MAPPING_DIR

indexes=(
    "Ath_Chame1"
    "Ath_Macpha1"
    "Ath_Parch1"
    "Ath_Phapo1"
    "Ath_Sorhu1"
    "Ath_Truan1"
    "Ath_Zalva1"
)
Mocks=(
    "$READS_DIR/4426_H_run639_CAGGTGTC_S134_L006_R1_001.fastq.trimmed.fastq"
    "$READS_DIR/4426_P_run639_AGATGCTA_S143_L007_R1_001.fastq.trimmed.fastq"
    "$READS_DIR/4426_X_run639_GCCGAGCT_S153_L008_R1_001.fastq.trimmed.fastq"
)
Chame1=(
    "$READS_DIR/4426_A_run639_TTCGACTA_S125_L005_R1_001.fastq.trimmed.fastq"
    "$READS_DIR/4426_I_run639_GCCGCATG_S135_L006_R1_001.fastq.trimmed.fastq"
    "$READS_DIR/4426_Q_run639_AGACACCT_S144_L007_R1_001.fastq.trimmed.fastq"
)
Macpha1=(
    "$READS_DIR/4426_E_run639_TTACCGAA_S129_L005_R1_001.fastq.trimmed.fastq"
    "$READS_DIR/4426_M_run639_CCGACGTT_S130_L005_R1_001.fastq.trimmed.fastq"
    "$READS_DIR/4426_U_run639_ACACAGTG_S150_L008_R1_001.fastq.trimmed.fastq"
)
Parch1=(
    "$READS_DIR/4426_B_run639_GTCTTACA_S126_L005_R1_001.fastq.trimmed.fastq"
    "$READS_DIR/4426_J_run639_CTGCACAT_S136_L006_R1_001.fastq.trimmed.fastq"
    "$READS_DIR/4426_R_run639_TTATAGCC_S145_L007_R1_001.fastq.trimmed.fastq"
)
Phapo1=(
    "$READS_DIR/4426_C_run639_GCATGAAC_S127_L005_R1_001.fastq.trimmed.fastq"
    "$READS_DIR/4426_K_run639_CGTGTCTA_S137_L006_R1_001.fastq.trimmed.fastq"
    "$READS_DIR/4426_S_run639_GTCTATAA_S146_L007_R1_001.fastq.trimmed.fastq"
)
Sorhu1=(
    "$READS_DIR/4426_F_run639_GTCTTACA_S139_L006_R1_001.fastq.trimmed.fastq"
    "$READS_DIR/4426_V_run639_TGTTCGAG_S151_L008_R1_001.fastq.trimmed.fastq"
    "$READS_DIR/4426_N_run639_CAATGAGC_S141_L007_R1_001.fastq.trimmed.fastq"
)
Truan1=(
    "$READS_DIR/4426_D_run639_CTCAATGA_S128_L005_R1_001.fastq.trimmed.fastq"
    "$READS_DIR/4426_L_run639_ATACGGAG_S138_L006_R1_001.fastq.trimmed.fastq"
    "$READS_DIR/4426_T_run639_GATCTTGC_S154_L008_R1_001.fastq.trimmed.fastq"
)
Zalva1=(
    "$READS_DIR/4426_G_run639_CGATCACA_S133_L006_R1_001.fastq.trimmed.fastq"
    "$READS_DIR/4426_O_run639_ATACGGAG_S142_L007_R1_001.fastq.trimmed.fastq"
    "$READS_DIR/4426_W_run639_GTCTATAA_S152_L008_R1_001.fastq.trimmed.fastq"
)

for mock in "${Mocks[@]}"; do
    for index in "${indexes[@]}"; do
        letter=$(python -c "print('$mock'.split('/')[-1].split('_')[1])")
        hisat2 -p $CORENUM -x $INDEX_DIR/$index -U $mock -S $MAPPING_DIR/Mock_$letter\_$index.sam
        samtools sort -@ $CORENUM -o $MAPPING_DIR/Mock_$letter\_$index.bam $MAPPING_DIR/Mock_$letter\_$index.sam
	samtools index $MAPPING_DIR/Mock_$letter\_$index.bam
	rm $MAPPING_DIR/Mock_$letter\_$index.sam
	cat $FILES_DIR/$(python -c "print('$index'.split('_')[0])").gtf $FILES_DIR/$(python -c "print('$index'.split('_')[1])").gtf > $FILES_DIR/$index.gtf
    done
done

Fungi=(
    "Chame1"
    "Macpha1"
    "Parch1"
    "Phapo1"
    "Sorhu1"
    "Truan1"
    "Zalva1"
)



for fung in "${Fungi[@]}"; do
    index=Ath_$fung
    arr=$fung[@]
    for sample in "${!arr}"; do
        letter=$(python -c "print('$sample'.split('/')[-1].split('_')[1])")
        hisat2 -p $CORENUM -x $INDEX_DIR/$index -U $sample -S $MAPPING_DIR/$fung\_$letter\_$index.sam
        samtools sort -@ $CORENUM -o $MAPPING_DIR/$fung\_$letter\_$index.bam $MAPPING_DIR/$fung\_$letter\_$index.sam
        samtools index $MAPPING_DIR/$fung\_$letter\_$index.bam
    done
done


