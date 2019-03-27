# contents of the directory at the start

```
.
|-- MAFBRA00707_R1.fastq
|-- MAFBRA00707_R2.fastq
|-- MAF.sh
|-- Mcanettii_R1.fastq.gz
|-- Mcanettii_R2.fastq.gz
|-- NC000962_3.fasta
|-- NC000962_3.gbk
`-- NexteraPE-PE.fa

```
# gzip these files


```


gzip -dc MAFBRA00707_R1.fastq.gz > MAFBRA00707_R1.fastq


gzip -dc MAFBRA00707_R2.fastq.gz > MAFBRA00707_R2.fastq

```



# trimmomatic 


```

java -jar /opt/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 MAFBRA00707_R1.fastq MAFBRA00707_R2.fastq MAFBRA00707_1_trimmed_paired.fastq MAFBRA00707_1_trimmed_unpaired.fastq MAFBRA00707_2_trimmed_paired.fastq MAFBRA00707_2_trimmed_unpaired.fastq ILLUMINACLIP:NexteraPE-PE.fa:2:30:10: LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36



```

## result of trimmomatic 

```
.
|-- after_trimmomatic.txt
|-- MAFBRA00707_1_trimmed_paired.fastq
|-- MAFBRA00707_1_trimmed_unpaired.fastq
|-- MAFBRA00707_2_trimmed_paired.fastq
|-- MAFBRA00707_2_trimmed_unpaired.fastq
|-- MAFBRA00707_R1.fastq
|-- MAFBRA00707_R1.fastq.gz
|-- MAFBRA00707_R2.fastq
|-- MAFBRA00707_R2.fastq.gz
|-- MAF.sh
|-- Mcanettii_R1.fastq.gz
|-- Mcanettii_R2.fastq.gz
|-- NC000962_3.fasta
|-- NC000962_3.gbk
`-- NexteraPE-PE.fa

0 directories, 15 files


```


# bwa_index_reference_genome 

```
bwa index NC000962_3.fasta

```

## after 

```

.
|-- after_bwa_index.txt
|-- after_trimmomatic.txt
|-- MAFBRA00707_1_trimmed_paired.fastq
|-- MAFBRA00707_1_trimmed_unpaired.fastq
|-- MAFBRA00707_2_trimmed_paired.fastq
|-- MAFBRA00707_2_trimmed_unpaired.fastq
|-- MAFBRA00707_R1.fastq
|-- MAFBRA00707_R1.fastq.gz
|-- MAFBRA00707_R2.fastq
|-- MAFBRA00707_R2.fastq.gz
|-- MAF.sh
|-- Mcanettii_R1.fastq.gz
|-- Mcanettii_R2.fastq.gz
|-- NC000962_3.fasta
|-- NC000962_3.fasta.amb
|-- NC000962_3.fasta.ann
|-- NC000962_3.fasta.bwt
|-- NC000962_3.fasta.pac
|-- NC000962_3.fasta.sa
|-- NC000962_3.gbk
`-- NexteraPE-PE.fa

0 directories, 21 files

```

# map_and_generate_sam_file 


```
bwa mem -R "@RG\tID:MAFBRA00707\tSM:MAFBRA00707\tPL:Illumina" -M NC000962_3.fasta MAFBRA00707_1_trimmed_paired.fastq MAFBRA00707_2_trimmed_paired.fastq > MAFBRA00707.sam

```

## after 

```
.
|-- after_bwa_index.txt
|-- after_bwa_mem.txt
|-- after_trimmomatic.txt
|-- MAFBRA00707_1_trimmed_paired.fastq
|-- MAFBRA00707_1_trimmed_unpaired.fastq
|-- MAFBRA00707_2_trimmed_paired.fastq
|-- MAFBRA00707_2_trimmed_unpaired.fastq
|-- MAFBRA00707_R1.fastq
|-- MAFBRA00707_R1.fastq.gz
|-- MAFBRA00707_R2.fastq
|-- MAFBRA00707_R2.fastq.gz
|-- MAFBRA00707.sam
|-- MAF.sh
|-- Mcanettii_R1.fastq.gz
|-- Mcanettii_R2.fastq.gz
|-- NC000962_3.fasta
|-- NC000962_3.fasta.amb
|-- NC000962_3.fasta.ann
|-- NC000962_3.fasta.bwt
|-- NC000962_3.fasta.pac
|-- NC000962_3.fasta.sa
|-- NC000962_3.gbk
`-- NexteraPE-PE.fa

0 directories, 23 files

```

# samtools_faidx_reference_genome 

```
samtools faidx NC000962_3.fasta
```

## after 

```
.
|-- after_bwa_index.txt
|-- after_bwa_mem.txt
|-- after_sam_faidx.txt
|-- after_trimmomatic.txt
|-- MAFBRA00707_1_trimmed_paired.fastq
|-- MAFBRA00707_1_trimmed_unpaired.fastq
|-- MAFBRA00707_2_trimmed_paired.fastq
|-- MAFBRA00707_2_trimmed_unpaired.fastq
|-- MAFBRA00707_R1.fastq
|-- MAFBRA00707_R1.fastq.gz
|-- MAFBRA00707_R2.fastq
|-- MAFBRA00707_R2.fastq.gz
|-- MAFBRA00707.sam
|-- MAF.sh
|-- Mcanettii_R1.fastq.gz
|-- Mcanettii_R2.fastq.gz
|-- NC000962_3.fasta
|-- NC000962_3.fasta.amb
|-- NC000962_3.fasta.ann
|-- NC000962_3.fasta.bwt
|-- NC000962_3.fasta.fai
|-- NC000962_3.fasta.pac
|-- NC000962_3.fasta.sa
|-- NC000962_3.gbk
`-- NexteraPE-PE.fa

0 directories, 25 files

```



# convert_sam_file_to_bam_file 

```
samtools view -bt NC000962_3.fasta.fai MAFBRA00707.sam > MAFBRA00707.bam
```


## after 

```
.
|-- after_bwa_index.txt
|-- after_bwa_mem.txt
|-- after_sam_faidx.txt
|-- after_sam_to_bam.txt
|-- after_trimmomatic.txt
|-- MAFBRA00707_1_trimmed_paired.fastq
|-- MAFBRA00707_1_trimmed_unpaired.fastq
|-- MAFBRA00707_2_trimmed_paired.fastq
|-- MAFBRA00707_2_trimmed_unpaired.fastq
|-- MAFBRA00707.bam
|-- MAFBRA00707_R1.fastq
|-- MAFBRA00707_R1.fastq.gz
|-- MAFBRA00707_R2.fastq
|-- MAFBRA00707_R2.fastq.gz
|-- MAFBRA00707.sam
|-- MAF.sh
|-- Mcanettii_R1.fastq.gz
|-- Mcanettii_R2.fastq.gz
|-- NC000962_3.fasta
|-- NC000962_3.fasta.amb
|-- NC000962_3.fasta.ann
|-- NC000962_3.fasta.bwt
|-- NC000962_3.fasta.fai
|-- NC000962_3.fasta.pac
|-- NC000962_3.fasta.sa
|-- NC000962_3.gbk
`-- NexteraPE-PE.fa

0 directories, 27 files

```


# sort_bam_file 


```
samtools sort MAFBRA00707.bam -o MAFBRA00707.sorted.bam
```

## after 

```
.
|-- after_bwa_index.txt
|-- after_bwa_mem.txt
|-- after_sam_faidx.txt
|-- after_sam_to_bam.txt
|-- after_samtools_sort.txt
|-- after_trimmomatic.txt
|-- MAFBRA00707_1_trimmed_paired.fastq
|-- MAFBRA00707_1_trimmed_unpaired.fastq
|-- MAFBRA00707_2_trimmed_paired.fastq
|-- MAFBRA00707_2_trimmed_unpaired.fastq
|-- MAFBRA00707.bam
|-- MAFBRA00707_R1.fastq
|-- MAFBRA00707_R1.fastq.gz
|-- MAFBRA00707_R2.fastq
|-- MAFBRA00707_R2.fastq.gz
|-- MAFBRA00707.sam
|-- MAFBRA00707.sorted.bam
|-- MAF.sh
|-- Mcanettii_R1.fastq.gz
|-- Mcanettii_R2.fastq.gz
|-- NC000962_3.fasta
|-- NC000962_3.fasta.amb
|-- NC000962_3.fasta.ann
|-- NC000962_3.fasta.bwt
|-- NC000962_3.fasta.fai
|-- NC000962_3.fasta.pac
|-- NC000962_3.fasta.sa
|-- NC000962_3.gbk
`-- NexteraPE-PE.fa

0 directories, 29 files

```



# samtools_index_sorted_bam 

```
samtools index MAFBRA00707.sorted.bam
```

## after 

```
.
|-- after_bwa_index.txt
|-- after_bwa_mem.txt
|-- after_sam_faidx.txt
|-- after_sam_to_bam.txt
|-- after_samtools_sort.txt
|-- after_trimmomatic.txt
|-- MAFBRA00707_1_trimmed_paired.fastq
|-- MAFBRA00707_1_trimmed_unpaired.fastq
|-- MAFBRA00707_2_trimmed_paired.fastq
|-- MAFBRA00707_2_trimmed_unpaired.fastq
|-- MAFBRA00707.bam
|-- MAFBRA00707_R1.fastq
|-- MAFBRA00707_R1.fastq.gz
|-- MAFBRA00707_R2.fastq
|-- MAFBRA00707_R2.fastq.gz
|-- MAFBRA00707.sam
|-- MAFBRA00707.sorted.bam
|-- MAF.sh
|-- Mcanettii_R1.fastq.gz
|-- Mcanettii_R2.fastq.gz
|-- NC000962_3.fasta
|-- NC000962_3.fasta.amb
|-- NC000962_3.fasta.ann
|-- NC000962_3.fasta.bwt
|-- NC000962_3.fasta.fai
|-- NC000962_3.fasta.pac
|-- NC000962_3.fasta.sa
|-- NC000962_3.gbk
`-- NexteraPE-PE.fa

0 directories, 29 files

```



# mapping_statistics 

```
samtools flagstat MAFBRA00707.sorted.bam > MAFBRA00707_stats.txt
```


# samtools_mpileup 

```
samtools mpileup -Q 23 -d 2000 -C 50 -ugf NC000962_3.fasta MAFBRA00707.sorted.bam | bcftools call -O v -vm -o MAFBRA00707.raw.vcf
```

## after 

```
.
|-- after_bwa_index.txt
|-- after_bwa_mem.txt
|-- after_sam_faidx.txt
|-- after_sam_to_bam.txt
|-- after_samtools_mpileup.txt
|-- after_samtools_sort.txt
|-- after_trimmomatic.txt
|-- MAFBRA00707_1_trimmed_paired.fastq
|-- MAFBRA00707_1_trimmed_unpaired.fastq
|-- MAFBRA00707_2_trimmed_paired.fastq
|-- MAFBRA00707_2_trimmed_unpaired.fastq
|-- MAFBRA00707.bam
|-- MAFBRA00707_R1.fastq
|-- MAFBRA00707_R1.fastq.gz
|-- MAFBRA00707_R2.fastq
|-- MAFBRA00707_R2.fastq.gz
|-- MAFBRA00707.raw.vcf
|-- MAFBRA00707.sam
|-- MAFBRA00707.sorted.bam
|-- MAFBRA00707.sorted.bam.bai
|-- MAFBRA00707_stats.txt
|-- MAF.sh
|-- Mcanettii_R1.fastq.gz
|-- Mcanettii_R2.fastq.gz
|-- NC000962_3.fasta
|-- NC000962_3.fasta.amb
|-- NC000962_3.fasta.ann
|-- NC000962_3.fasta.bwt
|-- NC000962_3.fasta.fai
|-- NC000962_3.fasta.pac
|-- NC000962_3.fasta.sa
|-- NC000962_3.gbk
`-- NexteraPE-PE.fa

0 directories, 33 files

```



# vcfutils_filter 



```
vcfutils.pl varFilter -d 10 -D 2000 MAFBRA00707.raw.vcf > MAFBRA00707.filt.vcf
```



# bgzip_filt_file 

```
bgzip -c MAFBRA00707.filt.vcf > MAFBRA00707.filt.vcf.gz
```



# run_tabix 

```
tabix -p vcf MAFBRA00707.filt.vcf.gz
```



# snpEff 

```
java -Xmx4g -jar /opt/snpEff/snpEff.jar -no-downstream -no-upstream -v -c /opt/snpEff/snpEff.config NC000962_3 MAFBRA00707.filt.vcf > MAFBRA00707.ann.vcf.gz
```

# velveth_assembly 

```
velveth MAFBRA00707_41 41 -fastq -shortPaired  MAFBRA00707_1_trimmed_paired.fastq MAFBRA00707_1_trimmed_unpaired.fastq -fastq -short MAFBRA00707_2_trimmed_paired.fastq MAFBRA00707_2_trimmed_unpaired.fastq
```

## after 

```
.
|-- after_bwa_index.txt
|-- after_bwa_mem.txt
|-- after_sam_faidx.txt
|-- after_sam_to_bam.txt
|-- after_samtools_mpileup.txt
|-- after_samtools_sort.txt
|-- after_snpeff.txt
|-- after_trimmomatic.txt
|-- after_velveth_41.txt
|-- MAFBRA00707_1_trimmed_paired.fastq
|-- MAFBRA00707_1_trimmed_unpaired.fastq
|-- MAFBRA00707_2_trimmed_paired.fastq
|-- MAFBRA00707_2_trimmed_unpaired.fastq
|-- MAFBRA00707_41
|   |-- Log
|   |-- Roadmaps
|   `-- Sequences
|-- MAFBRA00707.ann.vcf.gz
|-- MAFBRA00707.bam
|-- MAFBRA00707.filt.vcf
|-- MAFBRA00707.filt.vcf.gz
|-- MAFBRA00707.filt.vcf.gz.tbi
|-- MAFBRA00707_R1.fastq
|-- MAFBRA00707_R1.fastq.gz
|-- MAFBRA00707_R2.fastq
|-- MAFBRA00707_R2.fastq.gz
|-- MAFBRA00707.raw.vcf
|-- MAFBRA00707.sam
|-- MAFBRA00707.sorted.bam
|-- MAFBRA00707.sorted.bam.bai
|-- MAFBRA00707_stats.txt
|-- MAF.sh
|-- Mcanettii_R1.fastq.gz
|-- Mcanettii_R2.fastq.gz
|-- NC000962_3.fasta
|-- NC000962_3.fasta.amb
|-- NC000962_3.fasta.ann
|-- NC000962_3.fasta.bwt
|-- NC000962_3.fasta.fai
|-- NC000962_3.fasta.pac
|-- NC000962_3.fasta.sa
|-- NC000962_3.gbk
|-- NexteraPE-PE.fa
|-- snpEff_genes.txt
`-- snpEff_summary.html

1 directory, 44 files

```


# velvetg_produce_graph 

```
velvetg MAFBRA00707_41 -exp_cov auto -cov_cutoff auto
```



## after 

```
.
|-- after_bwa_index.txt
|-- after_bwa_mem.txt
|-- after_sam_faidx.txt
|-- after_sam_to_bam.txt
|-- after_samtools_mpileup.txt
|-- after_samtools_sort.txt
|-- after_snpeff.txt
|-- after_trimmomatic.txt
|-- after_velvetg_41.txt
|-- after_velveth_41.txt
|-- MAFBRA00707_1_trimmed_paired.fastq
|-- MAFBRA00707_1_trimmed_unpaired.fastq
|-- MAFBRA00707_2_trimmed_paired.fastq
|-- MAFBRA00707_2_trimmed_unpaired.fastq
|-- MAFBRA00707_41
|   |-- contigs.fa
|   |-- Graph2
|   |-- LastGraph
|   |-- Log
|   |-- PreGraph
|   |-- Roadmaps
|   |-- Sequences
|   `-- stats.txt
|-- MAFBRA00707.ann.vcf.gz
|-- MAFBRA00707.bam
|-- MAFBRA00707.filt.vcf
|-- MAFBRA00707.filt.vcf.gz
|-- MAFBRA00707.filt.vcf.gz.tbi
|-- MAFBRA00707_R1.fastq
|-- MAFBRA00707_R1.fastq.gz
|-- MAFBRA00707_R2.fastq
|-- MAFBRA00707_R2.fastq.gz
|-- MAFBRA00707.raw.vcf
|-- MAFBRA00707.sam
|-- MAFBRA00707.sorted.bam
|-- MAFBRA00707.sorted.bam.bai
|-- MAFBRA00707_stats.txt
|-- MAF.sh
|-- Mcanettii_R1.fastq.gz
|-- Mcanettii_R2.fastq.gz
|-- NC000962_3.fasta
|-- NC000962_3.fasta.amb
|-- NC000962_3.fasta.ann
|-- NC000962_3.fasta.bwt
|-- NC000962_3.fasta.fai
|-- NC000962_3.fasta.pac
|-- NC000962_3.fasta.sa
|-- NC000962_3.gbk
|-- NexteraPE-PE.fa
|-- snpEff_genes.txt
`-- snpEff_summary.html

1 directory, 50 files
```

# assemblathon_stats 

```
assemblathon_stats.pl ./MAFBRA00707_41/contigs.fa > assemblathon_stats_41.txt
```



# velveth_assembly 

```
velveth MAFBRA00707_49 49 -fastq -shortPaired  MAFBRA00707_1_trimmed_paired.fastq MAFBRA00707_1_trimmed_unpaired.fastq -fastq -short MAFBRA00707_2_trimmed_paired.fastq MAFBRA00707_2_trimmed_unpaired.fastq
```



# velvetg_produce_graph 

```
velvetg MAFBRA00707_49 -exp_cov auto -cov_cutoff auto
```



# assemblathon_stats 

```
assemblathon_stats.pl ./MAFBRA00707_49/contigs.fa > assemblathon_stats_49.txt
```



# velveth_assembly 

```
velveth MAFBRA00707_55 55 -fastq -shortPaired  MAFBRA00707_1_trimmed_paired.fastq MAFBRA00707_1_trimmed_unpaired.fastq -fastq -short MAFBRA00707_2_trimmed_paired.fastq MAFBRA00707_2_trimmed_unpaired.fastq
```

# velvetg_produce_graph 

```
velvetg MAFBRA00707_55 -exp_cov auto -cov_cutoff auto
```


# assemblathon_stats 

```
assemblathon_stats.pl ./MAFBRA00707_55/contigs.fa > assemblathon_stats_55.txt
```



# comparative assemblathon stats 

```
assemblathon_stats.pl ./MAFBRA00707_41/contigs.fa

assemblathon_stats.pl ./MAFBRA00707_49/contigs.fa


assemblathon_stats.pl ./MAFBRA00707_55/contigs.fa  
```


# I ran the scala functions best_assemblathon_stats from https://github.com/BioDragao/tese/blob/master/analysis/Step4.sc#L850
# to find out the best assemblathon stats. Here's the result

# Highest quality k_mer : 55 


# abacas_align_contigs 

```
cd MAFBRA00707_55 &&  cp ../NC000962_3.fasta ./ && abacas.pl -r ../NC000962_3.fasta -q contigs.fa -p promer -b -d -a && cd ..
```

## after 

```
.
|-- after_abacas_align.txt
|-- after_bwa_index.txt
|-- after_bwa_mem.txt
|-- after_sam_faidx.txt
|-- after_sam_to_bam.txt
|-- after_samtools_mpileup.txt
|-- after_samtools_sort.txt
|-- after_snpeff.txt
|-- after_trimmomatic.txt
|-- after_velvetg_41.txt
|-- after_velvetg_49.txt
|-- after_velvetg_55.txt
|-- after_velveth_41.txt
|-- after_velveth_49.txt
|-- after_velveth_55.txt
|-- assemblathon_stats_41.txt
|-- assemblathon_stats_49.txt
|-- assemblathon_stats_55.txt
|-- MAFBRA00707_1_trimmed_paired.fastq
|-- MAFBRA00707_1_trimmed_unpaired.fastq
|-- MAFBRA00707_2_trimmed_paired.fastq
|-- MAFBRA00707_2_trimmed_unpaired.fastq
|-- MAFBRA00707_41
|   |-- contigs.fa
|   |-- Graph2
|   |-- LastGraph
|   |-- Log
|   |-- PreGraph
|   |-- Roadmaps
|   |-- Sequences
|   `-- stats.txt
|-- MAFBRA00707_49
|   |-- contigs.fa
|   |-- Graph2
|   |-- LastGraph
|   |-- Log
|   |-- PreGraph
|   |-- Roadmaps
|   |-- Sequences
|   `-- stats.txt
|-- MAFBRA00707_55
|   |-- contigs.fa
|   |-- contigs.fa_NC000962_3.fasta.bin
|   |-- contigs.fa_NC000962_3.fasta.contigsInbin.fas
|   |-- contigs.fa_NC000962_3.fasta.crunch
|   |-- contigs.fa_NC000962_3.fasta.fasta
|   |-- contigs.fa_NC000962_3.fasta.gaps
|   |-- contigs.fa_NC000962_3.fasta.gaps.stats
|   |-- contigs.fa_NC000962_3.fasta.gaps.tab
|   |-- contigs.fa_NC000962_3.fasta.tab
|   |-- Graph2
|   |-- LastGraph
|   |-- Log
|   |-- NC000962_3.fasta
|   |-- PreGraph
|   |-- Roadmaps
|   |-- Sequences
|   |-- stats.txt
|   `-- unused_contigs.out
|-- MAFBRA00707.ann.vcf.gz
|-- MAFBRA00707.bam
|-- MAFBRA00707.filt.vcf
|-- MAFBRA00707.filt.vcf.gz
|-- MAFBRA00707.filt.vcf.gz.tbi
|-- MAFBRA00707_R1.fastq
|-- MAFBRA00707_R1.fastq.gz
|-- MAFBRA00707_R2.fastq
|-- MAFBRA00707_R2.fastq.gz
|-- MAFBRA00707.raw.vcf
|-- MAFBRA00707.sam
|-- MAFBRA00707.sorted.bam
|-- MAFBRA00707.sorted.bam.bai
|-- MAFBRA00707_stats.txt
|-- MAF.sh
|-- Mcanettii_R1.fastq.gz
|-- Mcanettii_R2.fastq.gz
|-- NC000962_3.fasta
|-- NC000962_3.fasta.amb
|-- NC000962_3.fasta.ann
|-- NC000962_3.fasta.bwt
|-- NC000962_3.fasta.fai
|-- NC000962_3.fasta.pac
|-- NC000962_3.fasta.sa
|-- NC000962_3.gbk
|-- NexteraPE-PE.fa
|-- snpEff_genes.txt
`-- snpEff_summary.html


```

# prokka_annotation 

```
cd ./MAFBRA00707_55 && prokka --outdir ./MAFBRA00707_prokka --prefix MAFBRA00707 contigs.fa_NC000962_3.fasta.fasta && cd ..

```

## after 

```
.
|-- after_abacas_align.txt
|-- after_bwa_index.txt
|-- after_bwa_mem.txt
|-- after_prokka_annotation.txt
|-- after_sam_faidx.txt
|-- after_sam_to_bam.txt
|-- after_samtools_mpileup.txt
|-- after_samtools_sort.txt
|-- after_snpeff.txt
|-- after_trimmomatic.txt
|-- after_velvetg_41.txt
|-- after_velvetg_49.txt
|-- after_velvetg_55.txt
|-- after_velveth_41.txt
|-- after_velveth_49.txt
|-- after_velveth_55.txt
|-- assemblathon_stats_41.txt
|-- assemblathon_stats_49.txt
|-- assemblathon_stats_55.txt
|-- MAFBRA00707_1_trimmed_paired.fastq
|-- MAFBRA00707_1_trimmed_unpaired.fastq
|-- MAFBRA00707_2_trimmed_paired.fastq
|-- MAFBRA00707_2_trimmed_unpaired.fastq
|-- MAFBRA00707_41
|   |-- contigs.fa
|   |-- Graph2
|   |-- LastGraph
|   |-- Log
|   |-- PreGraph
|   |-- Roadmaps
|   |-- Sequences
|   `-- stats.txt
|-- MAFBRA00707_49
|   |-- contigs.fa
|   |-- Graph2
|   |-- LastGraph
|   |-- Log
|   |-- PreGraph
|   |-- Roadmaps
|   |-- Sequences
|   `-- stats.txt
|-- MAFBRA00707_55
|   |-- contigs.fa
|   |-- contigs.fa_NC000962_3.fasta.bin
|   |-- contigs.fa_NC000962_3.fasta.contigsInbin.fas
|   |-- contigs.fa_NC000962_3.fasta.crunch
|   |-- contigs.fa_NC000962_3.fasta.fasta
|   |-- contigs.fa_NC000962_3.fasta.gaps
|   |-- contigs.fa_NC000962_3.fasta.gaps.stats
|   |-- contigs.fa_NC000962_3.fasta.gaps.tab
|   |-- contigs.fa_NC000962_3.fasta.tab
|   |-- Graph2
|   |-- LastGraph
|   |-- Log
|   |-- MAFBRA00707_prokka
|   |   |-- errorsummary.val
|   |   |-- MAFBRA00707.ecn
|   |   |-- MAFBRA00707.err
|   |   |-- MAFBRA00707.faa
|   |   |-- MAFBRA00707.ffn
|   |   |-- MAFBRA00707.fixedproducts
|   |   |-- MAFBRA00707.fna
|   |   |-- MAFBRA00707.fsa
|   |   |-- MAFBRA00707.gbf
|   |   |-- MAFBRA00707.gff
|   |   |-- MAFBRA00707.log
|   |   |-- MAFBRA00707.sqn
|   |   |-- MAFBRA00707.tbl
|   |   |-- MAFBRA00707.tsv
|   |   |-- MAFBRA00707.txt
|   |   `-- MAFBRA00707.val
|   |-- NC000962_3.fasta
|   |-- PreGraph
|   |-- Roadmaps
|   |-- Sequences
|   |-- stats.txt
|   `-- unused_contigs.out
|-- MAFBRA00707.ann.vcf.gz
|-- MAFBRA00707.bam
|-- MAFBRA00707.filt.vcf
|-- MAFBRA00707.filt.vcf.gz
|-- MAFBRA00707.filt.vcf.gz.tbi
|-- MAFBRA00707_R1.fastq
|-- MAFBRA00707_R1.fastq.gz
|-- MAFBRA00707_R2.fastq
|-- MAFBRA00707_R2.fastq.gz
|-- MAFBRA00707.raw.vcf
|-- MAFBRA00707.sam
|-- MAFBRA00707.sorted.bam
|-- MAFBRA00707.sorted.bam.bai
|-- MAFBRA00707_stats.txt
|-- MAF.sh
|-- Mcanettii_R1.fastq.gz
|-- Mcanettii_R2.fastq.gz
|-- NC000962_3.fasta
|-- NC000962_3.fasta.amb
|-- NC000962_3.fasta.ann
|-- NC000962_3.fasta.bwt
|-- NC000962_3.fasta.fai
|-- NC000962_3.fasta.pac
|-- NC000962_3.fasta.sa
|-- NC000962_3.gbk
|-- NexteraPE-PE.fa
|-- _resources
|   |-- BRIG-0.95-dist
|   |   |-- bin
|   |   |-- BRIG.jar
|   |   |-- BRIGMANUAL.pdf
|   |   |-- cgview
|   |   |   |-- ca
|   |   |   |   `-- ualberta
|   |   |   |       `-- stothard
|   |   |   |           `-- cgview
|   |   |   |               |-- Cgview.class
|   |   |   |               |-- CgviewConstants.class
|   |   |   |               |-- CgviewFactory$1.class
|   |   |   |               |-- CgviewFactory$2.class
|   |   |   |               |-- CgviewFactory$3.class
|   |   |   |               |-- CgviewFactory$4.class
|   |   |   |               |-- CgviewFactory.class
|   |   |   |               |-- CgviewFactory$ElementDetails.class
|   |   |   |               |-- CgviewFactoryPtt.class
|   |   |   |               |-- CgviewFactoryTab.class
|   |   |   |               |-- CgviewHTMLDocument.class
|   |   |   |               |-- CgviewIO.class
|   |   |   |               |-- Feature.class
|   |   |   |               |-- FeatureRange.class
|   |   |   |               |-- FeatureSlot.class
|   |   |   |               |-- FeatureSlot$SortFeaturesByStart.class
|   |   |   |               |-- Feature$SortFeatureRangesByStart.class
|   |   |   |               |-- FileMover.class
|   |   |   |               |-- InnerLabel.class
|   |   |   |               |-- LabelBounds.class
|   |   |   |               |-- Label.class
|   |   |   |               |-- Legend.class
|   |   |   |               |-- LegendItem.class
|   |   |   |               |-- OuterLabel.class
|   |   |   |               |-- Plasmid.class
|   |   |   |               |-- SeriesImage.class
|   |   |   |               |-- SortLabelsByForceLabel.class
|   |   |   |               |-- SortLabelsByRadians.class
|   |   |   |               |-- SortLabelsByRadiansShift.class
|   |   |   |               |-- SortLabelsByRadius.class
|   |   |   |               `-- SortSeriesImageByZoomCenter.class
|   |   |   |-- cgview.jar
|   |   |   |-- cgview_xml_builder
|   |   |   |   |-- cgview_xml_builder.pl
|   |   |   |   |-- README
|   |   |   |   `-- test.sh
|   |   |   |-- includes
|   |   |   |   |-- as_png.png
|   |   |   |   |-- as_svg.png
|   |   |   |   |-- expand_in_g.png
|   |   |   |   |-- expand_in.png
|   |   |   |   |-- expand_out_g.png
|   |   |   |   |-- expand_out.png
|   |   |   |   |-- full_g.png
|   |   |   |   |-- full.png
|   |   |   |   |-- help.html
|   |   |   |   |-- help.png
|   |   |   |   |-- help_png.html
|   |   |   |   |-- info.js
|   |   |   |   |-- move_back_g.png
|   |   |   |   |-- move_back.png
|   |   |   |   |-- move_forward_g.png
|   |   |   |   |-- move_forward.png
|   |   |   |   |-- overlib.js
|   |   |   |   `-- stylesheet.css
|   |   |   |-- lib
|   |   |   |   |-- batik-awt-util.jar
|   |   |   |   |-- batik-dom.jar
|   |   |   |   |-- batik-svggen.jar
|   |   |   |   |-- batik-util.jar
|   |   |   |   |-- batik-xml.jar
|   |   |   |   |-- commons-lang-2.0.jar
|   |   |   |   |-- jargs.jar
|   |   |   |   `-- xercesImpl.jar
|   |   |   |-- manifestinfo
|   |   |   |-- README
|   |   |   `-- src
|   |   |       `-- ca
|   |   |           `-- ualberta
|   |   |               `-- stothard
|   |   |                   `-- cgview
|   |   |                       |-- Cgview.class
|   |   |                       |-- CgviewConstants.class
|   |   |                       |-- CgviewConstants.java
|   |   |                       |-- CgviewFactory$1.class
|   |   |                       |-- CgviewFactory$2.class
|   |   |                       |-- CgviewFactory$3.class
|   |   |                       |-- CgviewFactory$4.class
|   |   |                       |-- CgviewFactory.class
|   |   |                       |-- CgviewFactory$ElementDetails.class
|   |   |                       |-- CgviewFactory.java
|   |   |                       |-- CgviewFactoryPtt.class
|   |   |                       |-- CgviewFactoryPtt.java
|   |   |                       |-- CgviewFactoryTab.class
|   |   |                       |-- CgviewFactoryTab.java
|   |   |                       |-- CgviewHTMLDocument.class
|   |   |                       |-- CgviewHTMLDocument.java
|   |   |                       |-- CgviewIO.class
|   |   |                       |-- CgviewIO.java
|   |   |                       |-- Cgview.java
|   |   |                       |-- Feature.class
|   |   |                       |-- Feature.java
|   |   |                       |-- FeatureRange.class
|   |   |                       |-- FeatureRange.java
|   |   |                       |-- FeatureSlot.class
|   |   |                       |-- FeatureSlot.java
|   |   |                       |-- FeatureSlot$SortFeaturesByStart.class
|   |   |                       |-- Feature$SortFeatureRangesByStart.class
|   |   |                       |-- FileMover.class
|   |   |                       |-- FileMover.java
|   |   |                       |-- includes
|   |   |                       |   |-- as_png.png
|   |   |                       |   |-- as_svg.png
|   |   |                       |   |-- expand_in_g.png
|   |   |                       |   |-- expand_in.png
|   |   |                       |   |-- expand_out_g.png
|   |   |                       |   |-- expand_out.png
|   |   |                       |   |-- full_g.png
|   |   |                       |   |-- full.png
|   |   |                       |   |-- help.html
|   |   |                       |   |-- help.png
|   |   |                       |   |-- help_png.html
|   |   |                       |   |-- info.js
|   |   |                       |   |-- move_back_g.png
|   |   |                       |   |-- move_back.png
|   |   |                       |   |-- move_forward_g.png
|   |   |                       |   |-- move_forward.png
|   |   |                       |   |-- overlib.js
|   |   |                       |   `-- stylesheet.css
|   |   |                       |-- InnerLabel.class
|   |   |                       |-- InnerLabel.java
|   |   |                       |-- LabelBounds.class
|   |   |                       |-- LabelBounds.java
|   |   |                       |-- Label.class
|   |   |                       |-- Label.java
|   |   |                       |-- Legend.class
|   |   |                       |-- LegendItem.class
|   |   |                       |-- LegendItem.java
|   |   |                       |-- Legend.java
|   |   |                       |-- OuterLabel.class
|   |   |                       |-- OuterLabel.java
|   |   |                       |-- Plasmid.class
|   |   |                       |-- Plasmid.java
|   |   |                       |-- SeriesImage.class
|   |   |                       |-- SeriesImage.java
|   |   |                       |-- SortLabelsByForceLabel.class
|   |   |                       |-- SortLabelsByForceLabel.java
|   |   |                       |-- SortLabelsByRadians.class
|   |   |                       |-- SortLabelsByRadians.java
|   |   |                       |-- SortLabelsByRadiansShift.class
|   |   |                       |-- SortLabelsByRadiansShift.java
|   |   |                       |-- SortLabelsByRadius.class
|   |   |                       |-- SortLabelsByRadius.java
|   |   |                       |-- SortSeriesImageByZoomCenter.class
|   |   |                       `-- SortSeriesImageByZoomCenter.java
|   |   |-- COPYING.txt
|   |   |-- default-BRIG.xml
|   |   |-- errorlog.xml
|   |   |-- lib
|   |   |   |-- jdom.jar
|   |   |   `-- swing-layout-1.0.4.jar
|   |   |-- profiles
|   |   |   |-- Eight_rings.xml
|   |   |   |-- Eight_rings.xml.jpg
|   |   |   |-- Fifteen_rings.xml
|   |   |   |-- Fifteen_rings.xml.jpg
|   |   |   |-- Five_rings.xml
|   |   |   |-- Five_rings.xml.jpg
|   |   |   |-- Four_rings.xml
|   |   |   |-- Four_rings.xml.jpg
|   |   |   |-- Nine_rings.xml
|   |   |   |-- Nine_rings.xml.jpg
|   |   |   |-- One_ring.xml
|   |   |   |-- One_ring.xml.jpg
|   |   |   |-- Seven_rings.xml
|   |   |   |-- Seven_rings.xml.jpg
|   |   |   |-- Six_rings.xml
|   |   |   |-- Six_rings.xml.jpg
|   |   |   |-- Ten_rings.xml
|   |   |   |-- Ten_rings.xml.jpg
|   |   |   |-- Three_rings.xml
|   |   |   |-- Three_rings.xml.jpg
|   |   |   |-- Twenty_rings.xml
|   |   |   |-- Twenty_rings.xml.jpg
|   |   |   |-- Two_rings.xml
|   |   |   |-- Two_rings.xml.jpg
|   |   |   `-- Two_ring.xml
|   |   |-- proteins.txt
|   |   `-- README.txt
|   |-- BRIG-0.95-dist.zip
|   `-- ITGE2017_Documentation.pdf
|-- snpEff_genes.txt
`-- snpEff_summary.html

23 directories, 273 files

```

# snippy_command 

```
snippy --cpus 4 --outdir MAFBRA00707 --ref ./NC000962_3.gbk --R1 ./MAFBRA00707_R1.fastq.gz --R2 ./MAFBRA00707_R2.fastq.gz
```

# ALL DONE



TODO:

Now we added the analysis of RD-Analyzer technique to better understand the results 

https://github.com/BioDragao/RD-Analyzer


NOTE: This processs **fails** at some point
```

$ python2.7 RD-Analyzer-extended.py NC000962_3.fasta MAFBRA00707_R1.fastq 

[root@localhost RD-Analyzer]# python2.7 RD-Analyzer-extended.py NC000962_3.fasta MAFBRA00707_R1.fastq 
Files exist in the output directory with the output prefix specified, refuse to overwrite...
[root@localhost RD-Analyzer]# ls
MAFBRA00707_R1.fastq  RD-Analyzer.depth        RD-Analyzer.py
MAFBRA00707_R2.fastq  RD-Analyzer-extended.py  README.md
NC000962_3.fasta      RD-Analyzer-md.pdf       Reference
[root@localhost RD-Analyzer]# rm *.depth
rm: remove regular empty file ‘RD-Analyzer.depth’? y
[root@localhost RD-Analyzer]# python2.7 RD-Analyzer-extended.py NC000962_3.fasta MAFBRA00707_R1.fastq 
[bwa_index] Pack FASTA... 0.04 sec
[bwa_index] Construct BWT for the packed sequence...
[bwa_index] 1.79 seconds elapse.
[bwa_index] Update BWT... 0.02 sec
[bwa_index] Pack forward-only FASTA... 0.02 sec
[bwa_index] Construct SA from BWT and Occ... 0.76 sec
[main] Version: 0.7.12-r1039
[main] CMD: bwa index NC000962_3.fasta
[main] Real time: 2.984 sec; CPU: 2.646 sec
[M::bwa_idx_load_from_disk] read 0 ALT contigs
[M::process] read 66700 sequences (10000238 bp)...
[M::process] read 66736 sequences (10000161 bp)...
[M::mem_process_seqs] Processed 66700 reads in 5.094 CPU sec, 5.504 real sec
[M::process] read 66700 sequences (10000119 bp)...
[M::mem_process_seqs] Processed 66736 reads in 5.751 CPU sec, 6.327 real sec
[M::process] read 66712 sequences (10000279 bp)...
[M::mem_process_seqs] Processed 66700 reads in 6.106 CPU sec, 7.652 real sec
[M::process] read 66752 sequences (10000101 bp)...
[M::mem_process_seqs] Processed 66712 reads in 5.457 CPU sec, 6.028 real sec
[M::process] read 66678 sequences (10000041 bp)...
[M::mem_process_seqs] Processed 66752 reads in 10.181 CPU sec, 12.538 real sec
[M::process] read 66728 sequences (10000008 bp)...
[M::mem_process_seqs] Processed 66678 reads in 6.220 CPU sec, 7.012 real sec
[M::process] read 66740 sequences (10000286 bp)...
[M::mem_process_seqs] Processed 66728 reads in 7.151 CPU sec, 8.094 real sec
[M::process] read 66696 sequences (10000005 bp)...
[M::mem_process_seqs] Processed 66740 reads in 7.773 CPU sec, 9.177 real sec
[M::process] read 66734 sequences (10000154 bp)...
[M::mem_process_seqs] Processed 66696 reads in 6.174 CPU sec, 7.514 real sec
[M::process] read 66756 sequences (10000218 bp)...
[M::mem_process_seqs] Processed 66734 reads in 5.812 CPU sec, 6.854 real sec
[M::process] read 66698 sequences (10000060 bp)...
[M::mem_process_seqs] Processed 66756 reads in 5.878 CPU sec, 6.613 real sec
[M::process] read 66742 sequences (10000037 bp)...
[M::mem_process_seqs] Processed 66698 reads in 5.741 CPU sec, 6.624 real sec
[M::process] read 66742 sequences (10000122 bp)...
[M::mem_process_seqs] Processed 66742 reads in 6.146 CPU sec, 7.741 real sec
[M::process] read 66702 sequences (10000295 bp)...
[M::mem_process_seqs] Processed 66742 reads in 5.649 CPU sec, 6.093 real sec
[M::process] read 66718 sequences (10000151 bp)...
[M::mem_process_seqs] Processed 66702 reads in 5.079 CPU sec, 5.507 real sec
[M::process] read 66728 sequences (10000273 bp)...
[M::mem_process_seqs] Processed 66718 reads in 5.714 CPU sec, 6.656 real sec
[M::process] read 66742 sequences (10000218 bp)...
[M::mem_process_seqs] Processed 66728 reads in 5.871 CPU sec, 6.859 real sec
[M::process] read 66716 sequences (10000175 bp)...
[M::mem_process_seqs] Processed 66742 reads in 5.641 CPU sec, 6.214 real sec
[M::process] read 66680 sequences (10000275 bp)...
[M::mem_process_seqs] Processed 66716 reads in 5.715 CPU sec, 6.120 real sec
[M::process] read 66746 sequences (10000031 bp)...
[M::mem_process_seqs] Processed 66680 reads in 5.536 CPU sec, 6.262 real sec
[M::process] read 66732 sequences (10000196 bp)...
[M::mem_process_seqs] Processed 66746 reads in 5.870 CPU sec, 6.626 real sec
[M::process] read 66700 sequences (10000183 bp)...
[M::mem_process_seqs] Processed 66732 reads in 6.040 CPU sec, 7.260 real sec
[M::process] read 66742 sequences (10000248 bp)...
[M::mem_process_seqs] Processed 66700 reads in 7.593 CPU sec, 8.896 real sec
[M::process] read 66710 sequences (10000227 bp)...
^[[B[M::mem_process_seqs] Processed 66742 reads in 15.385 CPU sec, 20.562 real sec
[M::process] read 66726 sequences (10000164 bp)...
[M::mem_process_seqs] Processed 66710 reads in 8.121 CPU sec, 9.490 real sec
[M::process] read 66730 sequences (10000226 bp)...
[M::mem_process_seqs] Processed 66726 reads in 7.136 CPU sec, 8.119 real sec
[M::process] read 66702 sequences (10000187 bp)...
[M::mem_process_seqs] Processed 66730 reads in 6.321 CPU sec, 7.006 real sec
[M::process] read 66764 sequences (10000263 bp)...
[M::mem_process_seqs] Processed 66702 reads in 6.992 CPU sec, 7.729 real sec
[M::process] read 66748 sequences (10000258 bp)...
[M::mem_process_seqs] Processed 66764 reads in 7.984 CPU sec, 9.178 real sec
[M::process] read 66682 sequences (10000004 bp)...
[M::mem_process_seqs] Processed 66748 reads in 8.584 CPU sec, 10.135 real sec
[M::process] read 66744 sequences (10000271 bp)...
[M::mem_process_seqs] Processed 66682 reads in 6.057 CPU sec, 6.894 real sec
[M::process] read 66744 sequences (10000088 bp)...
[M::mem_process_seqs] Processed 66744 reads in 6.776 CPU sec, 7.526 real sec
[M::process] read 66680 sequences (10000266 bp)...
[M::mem_process_seqs] Processed 66744 reads in 7.178 CPU sec, 8.276 real sec
[M::process] read 66694 sequences (10000182 bp)...
[M::mem_process_seqs] Processed 66680 reads in 9.496 CPU sec, 12.065 real sec
[M::process] read 66662 sequences (10000158 bp)...
[M::mem_process_seqs] Processed 66694 reads in 6.344 CPU sec, 7.129 real sec
[M::process] read 66706 sequences (10000196 bp)...
[M::mem_process_seqs] Processed 66662 reads in 5.735 CPU sec, 6.225 real sec
[M::process] read 66702 sequences (10000017 bp)...
[M::mem_process_seqs] Processed 66706 reads in 5.818 CPU sec, 6.306 real sec
[M::process] read 66630 sequences (10000057 bp)...
[M::mem_process_seqs] Processed 66702 reads in 5.970 CPU sec, 6.502 real sec
[M::process] read 66706 sequences (10000175 bp)...
[M::mem_process_seqs] Processed 66630 reads in 5.727 CPU sec, 6.461 real sec
[M::process] read 66712 sequences (10000118 bp)...
[M::mem_process_seqs] Processed 66706 reads in 6.452 CPU sec, 8.452 real sec
[M::process] read 66654 sequences (10000257 bp)...
[M::mem_process_seqs] Processed 66712 reads in 6.032 CPU sec, 6.867 real sec
[M::process] read 66714 sequences (10000049 bp)...
[M::mem_process_seqs] Processed 66654 reads in 5.880 CPU sec, 6.384 real sec
[M::process] read 66702 sequences (10000185 bp)...
[M::mem_process_seqs] Processed 66714 reads in 6.318 CPU sec, 7.052 real sec
[M::process] read 66694 sequences (10000115 bp)...
[M::mem_process_seqs] Processed 66702 reads in 7.318 CPU sec, 8.164 real sec
[M::process] read 66716 sequences (10000250 bp)...
[M::mem_process_seqs] Processed 66694 reads in 8.899 CPU sec, 11.963 real sec
[M::process] read 66722 sequences (10000212 bp)...
[M::mem_process_seqs] Processed 66716 reads in 5.749 CPU sec, 6.918 real sec
[M::process] read 66702 sequences (10000037 bp)...
[M::mem_process_seqs] Processed 66722 reads in 5.769 CPU sec, 6.699 real sec
[M::process] read 66712 sequences (10000256 bp)...
[M::mem_process_seqs] Processed 66702 reads in 5.250 CPU sec, 5.965 real sec
[M::process] read 66712 sequences (10000099 bp)...
[M::mem_process_seqs] Processed 66712 reads in 8.741 CPU sec, 10.416 real sec
[M::process] read 66698 sequences (10000213 bp)...
[M::mem_process_seqs] Processed 66712 reads in 6.917 CPU sec, 7.583 real sec
[M::process] read 66634 sequences (10000119 bp)...
[M::mem_process_seqs] Processed 66698 reads in 5.146 CPU sec, 5.664 real sec
[M::process] read 66710 sequences (10000233 bp)...
[M::mem_process_seqs] Processed 66634 reads in 5.823 CPU sec, 6.519 real sec
[M::process] read 66696 sequences (10000120 bp)...
[M::mem_process_seqs] Processed 66710 reads in 7.582 CPU sec, 8.993 real sec
[M::process] read 66644 sequences (10000091 bp)...
[M::mem_process_seqs] Processed 66696 reads in 9.348 CPU sec, 11.182 real sec
[M::process] read 66726 sequences (10000284 bp)...
[M::mem_process_seqs] Processed 66644 reads in 7.928 CPU sec, 9.280 real sec
[M::process] read 66716 sequences (10000171 bp)...
[M::mem_process_seqs] Processed 66726 reads in 10.552 CPU sec, 13.564 real sec
[M::process] read 66684 sequences (10000059 bp)...
[M::mem_process_seqs] Processed 66716 reads in 6.965 CPU sec, 7.811 real sec
[M::process] read 66706 sequences (10000007 bp)...
[M::mem_process_seqs] Processed 66684 reads in 5.520 CPU sec, 6.587 real sec
[M::process] read 66712 sequences (10000278 bp)...
[M::mem_process_seqs] Processed 66706 reads in 5.631 CPU sec, 6.246 real sec
[M::process] read 66662 sequences (10000139 bp)...
[M::mem_process_seqs] Processed 66712 reads in 7.044 CPU sec, 9.362 real sec
[M::process] read 66724 sequences (10000222 bp)...
[M::mem_process_seqs] Processed 66662 reads in 9.945 CPU sec, 12.118 real sec
[M::process] read 66706 sequences (10000093 bp)...
[M::mem_process_seqs] Processed 66724 reads in 7.885 CPU sec, 9.254 real sec
[M::process] read 66656 sequences (10000150 bp)...
[M::mem_process_seqs] Processed 66706 reads in 5.229 CPU sec, 5.606 real sec
[M::process] read 66714 sequences (10000044 bp)...
[M::mem_process_seqs] Processed 66656 reads in 5.450 CPU sec, 5.865 real sec
[M::process] read 66730 sequences (10000189 bp)...
[M::mem_process_seqs] Processed 66714 reads in 5.328 CPU sec, 5.737 real sec
[M::process] read 66700 sequences (10000158 bp)...
[M::mem_process_seqs] Processed 66730 reads in 5.318 CPU sec, 5.960 real sec
[M::process] read 66744 sequences (10000157 bp)...
[M::mem_process_seqs] Processed 66700 reads in 5.430 CPU sec, 5.839 real sec
[M::process] read 3456 sequences (518170 bp)...
[M::mem_process_seqs] Processed 66744 reads in 5.145 CPU sec, 5.524 real sec
[M::mem_process_seqs] Processed 3456 reads in 0.285 CPU sec, 0.330 real sec
[main] Version: 0.7.12-r1039
[main] CMD: bwa mem -R @RG\tID:RD\tSM:RD\tLB:RD\tPL:Illumina NC000962_3.fasta MAFBRA00707_R1.fastq
[main] Real time: 532.861 sec; CPU: 457.690 sec
[bam_sort] Use -T PREFIX / -o FILE to specify temporary and final output files
Usage: samtools sort [options...] [in.bam]
Options:
  -l INT     Set compression level, from 0 (uncompressed) to 9 (best)
  -m INT     Set maximum memory per thread; suffix K/M/G recognized [768M]
  -n         Sort by read name
  -o FILE    Write final output to FILE rather than standard output
  -T PREFIX  Write temporary files to PREFIX.nnnn.bam
  -@, --threads INT
             Set number of sorting and compression threads [1]
      --input-fmt-option OPT[=VAL]
               Specify a single input file format option in the form
               of OPTION or OPTION=VALUE
  -O, --output-fmt FORMAT[,OPT[=VAL]]...
               Specify output format (SAM, BAM, CRAM)
      --output-fmt-option OPT[=VAL]
               Specify a single output file format option in the form
               of OPTION or OPTION=VALUE
      --reference FILE
               Reference sequence FASTA FILE [null]
[E::hts_open_format] fail to open file './RD-Analyzer.sort.bam'
samtools index: failed to open "./RD-Analyzer.sort.bam": No such file or directory
[E::hts_open_format] fail to open file './RD-Analyzer.sort.bam'
samtools depth: Could not open "./RD-Analyzer.sort.bam": No such file or directory
Traceback (most recent call last):
  File "RD-Analyzer-extended.py", line 241, in <module>
    info = t.readRef(reference)
  File "RD-Analyzer-extended.py", line 116, in readRef
    info[trace] = [0.09, 0.5, tmp[3]]
IndexError: list index out of range


```
