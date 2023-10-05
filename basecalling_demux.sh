#!/bin/bash

# command name: basecalling_demux.sh
# author: Luciano Kalabric Silva
# institution: Oswaldo Cruz Foundation, Goncalo Moniz Institute, Bahia, Brazil
# last update: 21 SET 2023
# objetive: Basecalling and demux MinIon data
# Syntax: ./basecalling_demux.sh <parameters>
# Link: https://timkahlke.github.io/LongRead_tutorials/BS_G.html

# Validating arguments
if [[ $# -ne 2 ]]; then
    echo "Illegal number of parameters"
    echo "Syntax: basecalling_demux.sh <BASECALLER: -g/-d> <RUN_ID> <MODEL: fast/hac/sup>"
    exit 0    
fi
BASECALLER=$1
RUN_ID=$2
MODEL=$3

# INPUT dir
INPUT_DIR="${HOME}/data/fast5/${RUN_ID}"

# OUTPUT dir
echo "Prepering folders to the analysis..."
OUTPUT_DIR="${HOME}/bioinfo-results/${RUN_ID}/"
[ ! -d "${OUTPUT_DIR}" ] && mkdir -vp ${OUTPUT_DIR}
BASECALLDIR="${OUTPUT_DIR}/BASECALL"
DEMUXDIR="${OUTPUT_DIR}/DEMUX"
DEMUXCATDIR="${OUTPUT_DIR}/DEMUX_CAT"

# Default Fast min_qscore=8; Hac min_qscore=9; Sup min_qscore=10
QSCORE=9
LENGTH=100

# Parâmetro de otimização minimap2, samtools, racon e kraken2
THREADS="$(lscpu | grep 'CPU(s):' | awk '{print $2}' | sed -n '1p')"

# Basecalling
case $BASECALLER in
    -g)
        # Parâmetros Guppy basecaller (ONT)
        CONFIG="dna_r9.4.1_450bps_${MODEL}.cfg" #dna_r9.4.1_450bps_fast.cfg dna_r9.4.1_450bps_hac.cfg dna_r9.4.1_450bps_sup.cfg
        # Parâmetros para otimização do Guppy basecaller (ONT) obtidos pelo benckmark utilizadando o LAPTOP-Yale
        case $MODEL in
          fast)
            GPUPERDEVICE=4		
            CHUNCKSIZE=1000		
            CHUNKPERRUNNER=50
          ;;
          hac)
            GPUPERDEVICE=12		
            CHUNCKSIZE=2000		
            CHUNKPERRUNNER=256
          ;;
          sup)
            GPUPERDEVICE=12		
            CHUNCKSIZE=1000		
            CHUNKPERRUNNER=256
          ;;
          *)
            GPUPERDEVICE=4		
            CHUNCKSIZE=1000		
            CHUNKPERRUNNER=50
          ;;
        esac
        # Cria a pasta BASECALLDIR e faz o basecalling
        [ ! -d $BASECALLDIR ] && mkdir -vp $BASECALLDIR
        fi
        echo -e "Running guppy_basecaller...\n"
        # Comando para guppy_basecaller usando GPU
        guppy_basecaller -r -i ${INPUT_DIR} -s "${BASECALLDIR}" -c ${CONFIG} -x auto --min_qscore ${QSCORE} --gpu_runners_per_device ${GPUPERDEVICE} --chunk_size ${CHUNCKSIZE} --chunks_per_runner ${CHUNKPERRUNNER} --verbose_logs
        
        # Demultiplex, adapter removal & sem headcrop 18 para uso do cutadapt
        # Parâmetros Guppy barcoder (ONT)
        ARRANGEMENTS="barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg"
        [ ! -d $DEMUXDIR ] && mkdir -vp $DEMUXDIR
        echo -e "Running guppy_barcoder em ${BASECALLDIR}...\n"
        guppy_barcoder -r -i "${BASECALLDIR}/pass" -s ${DEMUXDIR} --arrangements_files ${ARRANGEMENTS} --require_barcodes_both_ends  --detect_mid_strand_barcodes --trim_barcodes  
        # Renomeei a pasta contendo as reads unclassified para barcode00 para análise
        # [ -d "${DEMUXDIR}/unclassified" ] && mv "${DEMUXDIR}/unclassified" "${DEMUXDIR}/barcode00"
        # Concatena todos arquivos .fastq de cada barcode em um arquivo .fastq único
        [ ! -d ${DEMUXCATDIR} ] && mkdir -vp ${DEMUXCATDIR}
        for i in $(find ${DEMUXDIR} -mindepth 1 -type d -name "barcode*" -exec basename {} \; | sort); do
          [ -d "${DEMUXDIR}/${i}" ] && cat ${DEMUXDIR}/${i}/*.fastq > "${DEMUXCATDIR}/${i}.fastq"
          # Gera o arquivo de log
          # echo "${i} $(grep -c "runid" ${DEMUXDIR}/${i}/*.fastq)" >> ${DEMUXDIR}/passed_reads.log
          echo "${i} `echo $(cat ${DEMUXDIR}/${i}/*.fastq|wc -l)/4|bc`"  >> ${DEMUXDIR}/passed_reads.log
          # echo "${i} $(grep -c "runid" ${DEMUXCATDIR}/${i}.fastq)" >> ${DEMUXCATDIR}/passed_reads.log
          echo "${i} `echo $(cat ${DEMUXCATDIR}/${i}.fastq|wc -l)/4|bc`"  >> ${DEMUXCATDIR}/passed_reads.log
        done
        exit 0
        ;;
    -d)
        CONFIG="dna_r10.4.1_e8.2_400bps_${MODEL}.cfg" #dna_r10.4.1_e8.2_400bps_fast.cfg dna_r10.4.1_e8.2_400bps_hac.cfg dna_r10.4.1_e8.2_400bps_sup.cfg
        dorado download --model ${CONFIG}@v4.1.0
        dorado basecaller ${CONFIG}@v4.1.0 $INPUT_DIR > $OUTPUT_DIR/calls.bam
        exit 0
        ;;
    *)
        echo "Invalid parameter"
        ;;
esac
    

