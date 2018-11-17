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

```

$ python2.7 RD-Analyzer-extended.py NC000962_3.fasta MAFBRA00707_R1.fastq 


```
