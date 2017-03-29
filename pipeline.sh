#/bin/bash
#$1 is the id of the sample
export id=$1
mkdir -p temp_files/aligned/
mkdir -p temp_files/scans/
mkdir -p temp_files/coords/
mkdir -p temp_files/fastas/

echo aligning
#bowtie2 is available here: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
#Bowtie2Index is a genome index created with bowtie2-index
bowtie2 -X 2000 -x ~/genomes/human/hg19/Bowtie2Index -1 fastqs/$id.1.trimmed.paired.fastq.gz -2 fastqs/$id.2.trimmed.paired.fastq.gz -S temp_files/aligned/$id.sam.gz

echo converting sams to bams
#samtools is available here: http://www.htslib.org/
samtools view -bS temp_files/aligned/$id.sam.gz > temp_files/aligned/$id.bam

echo converting BAM of aligned reads to BED file
#bedtools is available here: http://bedtools.readthedocs.io/en/latest/
bedtools bamtobed -i temp_files/aligned/$id.bam | gzip -c > temp_files/coords/$id.bed.gz

echo removing mitochondrial reads, expand around cut sites
# cut sites are start when +, end when -
gzip -dc temp_files/coords/$id.bed.gz | grep -v -P '^chrM' | awk 'BEGIN{FS="\t";OFS="\t"};{if ($6=="+"){print $1, $2-50, $2+50, $4} else{print $1, $3-50, $3+50, $4}} > temp_files/coords/$id.regions.bed

echo sort bed file
sort -k1,1 -k2,2n temp_files/coords/$id.regions.bed > temp_files/coords/$id.regions.sorted.bed

echo merge overlapping BED entries
#note that this has the effect of merging nearby regions and also counting duplicate reads as a single read. 
bedtools merge -i temp_files/coords/$id.regions.sorted.bed | gzip -c > temp_files/coords/$id.unique.bed.gz

echo convert bed to coordinate file
gzip -dc temp_files/coords/$id.unique.bed.gz | awk '{print $1":"$2"-"$3}' | gzip -c  > temp_files/coords/$id.coord.gz

echo removing invalid coordinate entries
#tossInvalidTwoBitCoords.py is available here https://github.com/Carldeboer/GKM-PCA and removes invalid intervals (e.g. off chromosome)
#chrom_noMT.sizes is a chrom.sizes file without the chrM entry
tossInvalidTwoBitCoords.py -s ~/genomes/human/hg19/chrom_noMT.sizes -i temp_files/coords/$id.coord.gz -o temp_files/coords/$id.coord.fixed

echo to fasta
#twoBitToFa is available from Kent tools (http://hgdownload.cse.ucsc.edu/admin/exe/)
#allChrs.2bit is a 2bit genome file of the corresponding genome version (
twoBitToFa ~/genomes/human/hg19/allChrs.2bit temp_files/fastas/$id.fa -seqList=temp_files/coords/$id.coord.fixed
gzip temp_files/fastas/$id.fa

echo scan for k-mers
# AMUSED is available from https://github.com/Carldeboer/AMUSED
AMUSED -nsz -ds -bc -ns -s 8 -q temp_files/fastas/$id.fa.gz -o temp_files/scans/$id.scan
gzip temp_files/scans/$id.scan

echo get frequencies column
gzip -dc  temp_files/scans/$id.scan.gz | awk 'BEGIN{FS="\t"; OFS="\t"};{print $2,$6}' | gzip -c > temp_files/scans/$id.freq.gz


#once these k-mer frequency columns have been collected for each sample, you can input them individually into R using read.table
