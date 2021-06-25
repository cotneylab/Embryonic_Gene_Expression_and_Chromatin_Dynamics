#!/bin/bash
#SBATCH --job-name=HMkallisto
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=125G
#SBATCH --mail-user=wentworth@uchc.edu
#SBATCH -o %j.out
#SBATCH -e %j.err
source ~/.bashrc_miniconda3
conda activate bustools
cd /home/FCAM/ewentworth/cotney/rawdata/human/embryonic_craniofacial/multiome/cs20_01/
module load kallisto
kb count -i /home/FCAM/ewentworth/cotney/genome/hg38/hg38_kallisto_index.idx -g /home/FCAM/ewentworth/cotney/genome/hg38/t2g.txt -x 10XV3 -o /home/FCAM/ewentworth/cotney/rawdata/human/embryonic_craniofacial/multiome/cs20_01/rna-out -t 1 --overwrite /home/FCAM/ewentworth/cotney/rawdata/human/embryonic_craniofacial/multiome/cs20_01/JC20011_GT21-00475_CCCAGCTTCT-GTTTGGTGTC_S11_L002_R1_001.fastq.gz /home/FCAM/ewentworth/cotney/rawdata/human/embryonic_craniofacial/multiome/cs20_01/JC20011_GT21-00475_CCCAGCTTCT-GTTTGGTGTC_S11_L002_R2_001.fastq.gz
bustools whitelist rna-out/output.bus -o whitelist.txt
kb count -i /home/FCAM/ewentworth/cotney/genome/hg38/hg38_kallisto_index.idx -g /home/FCAM/ewentworth/cotney/genome/hg38/t2g.txt -x 10XV3 -w whitelist.txt -o /home/FCAM/ewentworth/cotney/rawdata/human/embryonic_craniofacial/multiome/cs20_01/rna-out -t 1 --overwrite /home/FCAM/ewentworth/cotney/rawdata/human/embryonic_craniofacial/multiome/cs20_01/JC20011_GT21-00475_CCCAGCTTCT-GTTTGGTGTC_S11_L002_R1_001.fastq.gz /home/FCAM/ewentworth/cotney/rawdata/human/embryonic_craniofacial/multiome/cs20_01/JC20011_GT21-00475_CCCAGCTTCT-GTTTGGTGTC_S11_L002_R2_001.fastq.gz
