import argparse
import sys
import os
from subprocess import call

parser = argparse.ArgumentParser(prog='main.py')
parser.add_argument('--threads', help='Number of threads name')
parser.add_argument('--reference', help='Reference genome')
args = parser.parse_args()

chipseq_files = ['H3K27ac_aso1_rep1',
'H3K27ac_aso1_rep2',
'H3K27ac_aso2_rep1',
'H3K27ac_aso2_rep2',
'H3K27ac_control_rep1',
'H3K27ac_control_rep2',
'H3K36me3_aso2_rep1',
'H3K36me3_aso2_rep2',
'H3K36me3_control_rep1',
'H3K36me3_control_rep2',
'H3K9me3_aso1_rep1',
'H3K9me3_aso1_rep2',
'H3K9me3_aso2_rep1',
'H3K9me3_aso2_rep2',
'H3K9me3_control_rep1',
'H3K9me3_control_rep2']

input_files = [name + "_input" for name in chipseq_files]

if not os.path.exists("trimmed"):
	os.mkdir("trimmed")

# Preprocessing
for file in (chipseq_files + input_files):
	if call(f"fastp --in1  {file}.fastq  --trim_poly_g  --adapter_fasta  adapters.fa --out1 trimmed/{file}.fastq --thread {args.threads}  --disable_quality_filtering --poly_g_min_len 5" , shell=True)!=0:
		sys.exit(f"fastp {file} failed")

	if call(f"fastqc  {file}.fastq" , shell=True)!=0:
		sys.exit(f"fastqc {file} failed")

	if call(f"fastqc trimmed/{file}.fastq" , shell=True)!=0:
		sys.exit(f"fastqc trimmed/{file} failed")

	if call(f"(bowtie2 -p {args.threads}  --local -x {args.reference}  -U  trimmed/{file}.fastq -S  {file}.sam) 2>{file}_bowtie2_report" , shell=True)!=0:
		sys.exit(f"bowtie2 {file} failed")

	if call(f"samtools view -@ {args.threads}  -bS {file}.sam > {file}.bam", shell=True)!=0:
		sys.exit(f"samtools view {file} failed")
		
	os.remove(f"{file}.sam")
	
	if call(f"samtools sort -@  {args.threads}   {file}.bam > {file}.sorted.bam", shell=True)!=0:
		sys.exit(f"samtools sort {file} failed")
		
	if call(f"(samtools markdup -@  {args.threads}  -r -s  {file}.sorted.bam {file}.dedup.bam) 2> {file}_dedup_report", shell=True)!=0:
		sys.exit(f"samtools markdup {file} failed")
		
	if call(f"bedtools bamtobed -i {file}.dedup.bam > {file}.bed", shell=True)!=0:
		sys.exit(f"bedtools bamtobed  {file} failed")		

	if call(f"bedtools genomecov -ibam   {file}.dedup.bam  -bg >  {file}.bedgraph", shell=True)!=0:
		sys.exit(f"bedtools genomecov  {file} failed")
		
	if call(f"sort -k1,1 -k2,2n {file}.bedgraph >  {file}.sorted.bedgraph", shell=True)!=0:
		sys.exit(f"sort bedgraph {file} failed")
		
	if call(f"bedGraphToBigWig {file}.sorted.bedgraph  hg38.chrom.sizes  {file}.bw", shell=True)!=0:
		sys.exit(f"bedGraphToBigWig {file} failed")

# Sicer and Macs peak calling
for sig,inp in zip(chipseq_files,input_files):
	if call(f"sicer  -t {sig}.dedup.bam -c {inp}.dedup.bam -s hg38 --cpu {args.threads} --false_discovery_rate 0.05 > {sig}_sicer_report", shell=True)!=0:
		sys.exit(f"sicer {sig} failed")
		
	if call(f"sed 's/[\t]*$//' {sig}.dedup-W200-G600-FDR0.05-island.bed > {sig}.dedup-W200-G600-FDR0.05-island_corrected.bed", shell=True)!=0:
		sys.exit(f"sed {sig} failed")

	if call(f"bedtools intersect -v -a {sig}.dedup-W200-G600-FDR0.05-island_corrected.bed  -b hg38.blacklist.bed >  {sig}.sicer.bed", shell=True)!=0:
		sys.exit(f"bedtools intersect blacklist {sig} failed")
		
	if call(f"(macs2 callpeak -t {sig}.dedup.bam -c {inp}.dedup.bam -n {sig} --outdir {sig}  -f BAM --keep-dup all -g 3.0e9 --broad  -p 0.0005) 2> {sig}_macs_report", shell=True)!=0:
		sys.exit(f"macs2 callpeak {sig} failed")
		
	if call(f"bedtools intersect -v -a {sig}/{sig}_peaks.broadPeak -b hg38.blacklist.bed >  {sig}.macs.bed", shell=True)!=0:
		sys.exit(f"bedtools intersect blacklist {sig} macs failed")
		
# Homer peak calling
for file in (chipseq_files + input_files):
	if call(f"makeTagDirectory {file} {file}.dedup.bam -format sam" , shell=True)!=0:
		sys.exit(f"makeTagDirectory {file} failed")
		
for sig,inp in zip(chipseq_files,input_files):

	if call(f"findPeaks {sig} -style histone -i {inp} -fdr 0.05 -o {sig}_homer" , shell=True)!=0:
		sys.exit(f"findPeaks {sig} failed")
		
	if call(f"grep -v '^#' {sig}_homer > {sig}_homer.tmp" , shell=True)!=0:
		sys.exit(f"grep {sig} failed")	
	
	if call(f"perl pos2bedmod.pl {sig}_homer.tmp > {sig}_homer.bed" , shell=True)!=0:
		sys.exit(f"pos2bedmod {sig} failed")
		
	if call(f"bedtools intersect -v -a {sig}_homer.bed -b hg38.blacklist.bed >  {sig}.homer.bed", shell=True)!=0:
		sys.exit(f"bedtools intersect blacklist {sig} homer failed")
		




