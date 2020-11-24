#! /bin/bash

# Claire Vandiedonck - RNASeq M2 MEG - 2018
# Wednesday 5th December 2018
# Practical sessions 4 & 5:
# 	4. Mapping quality check and visualizing mapped reads
#	5. Counting reads per feature


#############
# Session 4
#############

#-------------------------------------------------------------------------
# 4. QC on bam files
#-------------------------------------------------------------------------

# Qualimap is a tool allowing to asses the quality of the mapped data, including coverage estimation and many other parameters
# It reads sorted.bam files and generates a folder containing a report on .html format
# Run it on Hypoxia and Normoxia data with the commands below where:
# 		bamqc means it will perform a qc on bam files
#		-bam sepcifies the format of input files, followed by the name of the input files
#		-c displays chromosome limits on graphics
#		--feature-file to use a feature annotation file, here the .gff file
# 		-outdir is used to specify the name of the output folder
qualimap bamqc -bam Normoxia_1_bowtie_mapping.sorted.bam -c --feature-file C_parapsilosis_ORFs.gff -outdir normoxia_qualimap
qualimap bamqc -bam Hypoxia_1_bowtie_mapping.sorted.bam -c --feature-file C_parapsilosis_ORFs.gff -outdir hypoxia_qualimap

# to look at the results, either double-click on the qualimapReport.html in the generated folder or use the command line below that opens the html file
firefox normoxia_qualimap/qualimapReport.html &
firefox hypoxia_qualimap/qualimapReport.html &
# Warning: the display of the letters and special caracters in the html file is by default "occidental".
# To have a proper display of all characters, open the menu on the top right hand corner.
# Click on "customize" (or "Personnaliser" in French) and select the text encoding icon.
# Slide it to the menu on the right. It now appears in your menu bar (or after clicking on the double arrow >>).
# Click on it and select "Unicode" instead of "occidental".

 # Interpret the outputs:
# 1) What are the proportions of reads properly mapped to the reference genome?
# 2) What do you think about the depth?
# 3) How do you explain that some reads were not mapped?
# 4) Is there a link between the fastqc read qualities and the bam qualities?

## ANSWERS:
## 1) Hypoxia: 91.28% of the reads are mapped, hence 8.72% unmapped.
## Normoxia: 89.69% of the reads are mapped, hence 10.31% unmapped (about 1 millions of reads)
## The duplication rate in both cases is slightly above 50%.
## In genic regions, the % of mapped reads and their number is similar in both conditions: 68% in the target region covered by 7.9 millions of reads (78%) for Hypoxia, versus 68% for Normoxia with 7.2 millions of reads (71%).
## 2) The depth of coverage in annotated regions is 36.47X with a standard deviation of 258.25 for Hypoxia, versus 31.79X with a sd of 111.25 for Normoxia.
## These depths are fully satisfaying for RNASeq. They are of the same range for the 8 chromosomes of CParapislosis except for the mitochondrial DNA
## better covered for Normoxia.
## 3) The unmapped reads may result from:
## - libraries issues (there are duplicates and we saw with the fastaqc that some k-mers were overrepresented)
## - amplification biases (~50% of duplicates)
## - sequencing errors
## - mapping errors that are only 1.35% for Hypoxia whereas 7.66% for Normoxia where more uncalled N bases were observed 
## - and mostly because reads may fall outside the known annotated regions => in the paper, they wanted to discover new features
## 4) With the fastqc reports, we saw that quality was really good for the Hypoxia sample whereas the quality was lower for Normoxia with two populations of reads, those with good Q scores and those, more numerous, with bad Q scores. These QC observations are in accordance with a higher error rate in the mapping of the Normoxia samples but they do not affect the depth of coverage in this dataset.

#-------------------------------------------------------------------------
# 5. Visualisation using IGV
#-------------------------------------------------------------------------

# IGV is a viewer for genomic data, including bam files once they are sorted and indexed
# Run IGV by typing the following command
igv
# load the genome in fasta format using the menu Genome/Load Genome from File/ 
# load the annotations as .GFF files using the menu File/Load from File
# load the sorted bam (sorted.bam.bai must be in the same directory as your sorted bam) also using the menu File/Load from File
# To visualize reads, you need to zoom a lot. You may do so by using the + in the top left hand corner or by sliding a gate on the genome with your mouse.
# You may also look at some specific genomic regions or genes by entering their coordinates or names in the navigator search tool.

# Is your RNASeq strand specific?
# Do you see any difference between hypoxia and normoxia files?
# Do your data fit the results of the paper?

## ANSWERS:
## Hypoxia: strand-sepcific with reads in the reverse complement orientation of the transcript, thus meaning dUTPs were incorporated for the second strand syntehsis
## Normoxia: not strand-sepcific, sequenced earlier with an older chemistry
## For each sample, we can detect new transcripts (ex: contig005807_C_parapislosis_CDC317 at 1668 kb in normoxia). It is impossible to asses whether the transcript is on the forward or reverse strand for Normoxia.
## In accordance with Qualimamp, there is on overall a better coverage for Normoxia than for Hypoxia...but these differences may also be due to differential expression between these conditions!


#############
# Session 5
#############

#-------------------------------------------------------------------------
# 6. Read counts per gene using BEDTOOLS on Unix or/and htseq count on Galaxie (biomina server)
#-------------------------------------------------------------------------

# Quantify your data by annotated features using BEDTOOLS with the following command
# multicov is to count the number of reads to multiple features
# -bams is to specify you are using bam files rather than sam files
# -bed is the option to specify the name of the annotation file, here the .gff file
# > is to direct the output to the file that you name
# (approximate time = 1 minute each)
bedtools multicov -bams Normoxia_1_bowtie_mapping.sorted.bam -bed C_parapsilosis_ORFs.gff > Normoxia_1_gene_counts.gff
bedtools multicov -bams Hypoxia_1_bowtie_mapping.sorted.bam -bed C_parapsilosis_ORFs.gff > Hypoxia_1_gene_counts.gff

# you may open them with less -SN as you used yesterday for the fastq data

# sed = stream editor, is a powerful tool in Unix to handle and edit text files
# here, you use it to delete all the characters from the beginning of each line up to ID= included, in order to only keep the last two columns, one with the gene name, the other with the read counts 
sed 's/^.*ID=//' Normoxia_1_gene_counts.gff > Normoxia_1_counts.tab
sed 's/^.*ID=//' Hypoxia_1_gene_counts.gff > Hypoxia_1_counts.tab

# What are the problems associated to this way of counting reads to features?
# Which other methods could have been used?
# Does it make sens to compare the two samples on this basis?

## ANSWERS:
## 1st problem: the counting mode for overlapping reads in ORFs. Which one should be chosen? Should we specify the strand? Here, we counted the read as long as at least one of its base did overlap the ORF which accounting for the strandedness. We could also use HTseq to specify different modes like uniion, strict intersection or empty intersection.
## 2nd problem: normalisation. The method we used doesn't account for the library size or for the numbers of mapped reads (and should we remove duplicates prior to the counting?). There were more reads mapped for Normoxia than for Hypoxia. The comparision between the two smaples should include a normalisation step before to adjust for the library size? We could also normalise according to the length of the genes: here, since we will want to compare the expression of each gene between two conditions, the size of the gene does not change between the two samples so this issue is less critical. Should we want to compare the expression between genes, that would be critical. It would also have an impact of we would like to look at co-expression networks.
##Students can compute the RPKM if they wish to but in any case with only one sample per condition the comparison would rely on statistics. For statistics, we need replicates!!!


