
```
# gzip these files

gzip -dc MAFBRA00707_R1.fastq.gz > MAFBRA00707_R1.fastq


gzip -dc MAFBRA00707_R2.fastq.gz > MAFBRA00707_R2.fastq



# trimmomatic <<<<<



java -jar /opt/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 MAFBRA00707_R1.fastq MAFBRA00707_R2.fastq MAFBRA00707_1_trimmed_paired.fastq MAFBRA00707_1_trimmed_unpaired.fastq MAFBRA00707_2_trimmed_paired.fastq MAFBRA00707_2_trimmed_unpaired.fastq ILLUMINACLIP:NexteraPE-PE.fa:2:30:10: LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36


# bwa_index_reference_genome <<<<<


bwa index NC000962_3.fasta



# map_and_generate_sam_file <<<<<


bwa mem -R "@RG\tID:MAFBRA00707\tSM:MAFBRA00707\tPL:Illumina" -M NC000962_3.fasta MAFBRA00707_1_trimmed_paired.fastq MAFBRA00707_2_trimmed_paired.fastq > MAFBRA00707.sam


# samtools_faidx_reference_genome <<<<<


samtools faidx NC000962_3.fasta




# convert_sam_file_to_bam_file <<<<<


samtools view -bt NC000962_3.fasta.fai MAFBRA00707.sam > MAFBRA00707.bam




# sort_bam_file <<<<<


samtools sort MAFBRA00707.bam -o MAFBRA00707.sorted.bam




# samtools_index_sorted_bam <<<<<


samtools index MAFBRA00707.sorted.bam




# mapping_statistics <<<<<


samtools flagstat MAFBRA00707.sorted.bam > MAFBRA00707_stats.txt




# samtools_mpileup <<<<<


samtools mpileup -Q 23 -d 2000 -C 50 -ugf NC000962_3.fasta MAFBRA00707.sorted.bam | bcftools call -O v -vm -o MAFBRA00707.raw.vcf




# vcfutils_filter<<<<<


vcfutils.pl varFilter -d 10 -D 2000 MAFBRA00707.raw.vcf > MAFBRA00707.filt.vcf




# bgzip_filt_file <<<<<


bgzip -c MAFBRA00707.filt.vcf > MAFBRA00707.filt.vcf.gz




# run_tabix <<<<<


tabix -p vcf MAFBRA00707.filt.vcf.gz




# snpEff <<<<<


java -Xmx4g -jar /opt/snpEff/snpEff.jar -no-downstream -no-upstream -v -c /opt/snpEff/snpEff.config NC000962_3 MAFBRA00707.filt.vcf > MAFBRA00707.ann.vcf.gz


# velveth_assembly <<<<<


velveth MAFBRA00707_41 41 -fastq -shortPaired  MAFBRA00707_1_trimmed_paired.fastq MAFBRA00707_1_trimmed_unpaired.fastq -fastq -short MAFBRA00707_2_trimmed_paired.fastq MAFBRA00707_2_trimmed_unpaired.fastq



# velvetg_produce_graph <<<<<


velvetg MAFBRA00707_41 -exp_cov auto -cov_cutoff auto




# assemblathon_stats <<<<<


assemblathon_stats.pl ./MAFBRA00707_41/contigs.fa > assemblathon_stats_41.txt





# velveth_assembly <<<<<


velveth MAFBRA00707_49 49 -fastq -shortPaired  MAFBRA00707_1_trimmed_paired.fastq MAFBRA00707_1_trimmed_unpaired.fastq -fastq -short MAFBRA00707_2_trimmed_paired.fastq MAFBRA00707_2_trimmed_unpaired.fastq




# velvetg_produce_graph <<<<<


velvetg MAFBRA00707_49 -exp_cov auto -cov_cutoff auto




# assemblathon_stats <<<<<


assemblathon_stats.pl ./MAFBRA00707_49/contigs.fa > assemblathon_stats_49.txt




# velveth_assembly <<<<<


velveth MAFBRA00707_55 55 -fastq -shortPaired  MAFBRA00707_1_trimmed_paired.fastq MAFBRA00707_1_trimmed_unpaired.fastq -fastq -short MAFBRA00707_2_trimmed_paired.fastq MAFBRA00707_2_trimmed_unpaired.fastq

# velvetg_produce_graph <<<<<


velvetg MAFBRA00707_55 -exp_cov auto -cov_cutoff auto



# assemblathon_stats <<<<<


assemblathon_stats.pl ./MAFBRA00707_55/contigs.fa > assemblathon_stats_55.txt




# comparative assemblathon stats 

assemblathon_stats.pl ./MAFBRA00707_41/contigs.fa


assemblathon_stats.pl ./MAFBRA00707_49/contigs.fa


assemblathon_stats.pl ./MAFBRA00707_55/contigs.fa  

# I ran the scala functions best_assemblathon_stats from https://github.com/BioDragao/tese/blob/master/analysis/Step4.sc#L850
# to find out the best assemblathon stats. Here's the result
# Highest quality k_mer : 55 


# abacas_align_contigs <<<<<


cd MAFBRA00707_55 &&  cp ../NC000962_3.fasta ./ && abacas.pl -r ../NC000962_3.fasta -q contigs.fa -p promer -b -d -a && cd ..


# prokka_annotation <<<<<


cd ./MAFBRA00707_55 && prokka --outdir ./MAFBRA00707_prokka --prefix MAFBRA00707 contigs.fa_NC000962_3.fasta.fasta && cd ..




# snippy_command <<<<<


snippy --cpus 4 --outdir MAFBRA00707 --ref ./NC000962_3.gbk --R1 ./MAFBRA00707_R1.fastq.gz --R2 ./MAFBRA00707_R2.fastq.gz


# ALL DONE

```
