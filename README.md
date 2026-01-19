Introduction
The objective of this assignment is to assemble a consensus genome of Salmonella enterica from Oxford Nanopore R10 sequencing reads, align the assembly to a reference genome obtained from NCBI, and finally, identify and visualize genomic differences. This workflow highlights both the strengths and limitations of long-read genome assembly and some genome comparison approaches.
Genome assembly is a core task in bioinformatics that aims to read and reconstruct an organism’s whole genome sequence from fragmented sequencing reads. Recent long-read sequencing technologies, such as Oxford Nanopore Technologies, improve the feasibility of de novo genome assembly for bacterial genomes (Salmonella enterica) [1]. Additionally, genome assembly remains challenging due to sequencing errors, repetitive genomic regions, costs and structural complexity.
Long-read sequencing offers a major advantage over short-read technologies by producing reads that can cover repetitive regions and complex genomic structures, resulting in better assembly contiguity [1]. However, long-read data could be more error-prone than short-read data, such as the insertion–deletion error. This affects consensus accuracy and complicates downstream analyses if not properly corrected during assembly and polishing steps [1].
Therefore, Flye assemblers have been developed to address the challenges of assembling long, error-prone reads by using repeat graphs for complex genomic structures [2]. For relatively small bacterial genomes like Salmonella enterica (5 Mb), only long-read assembly is often sufficient to produce highly contiguous assemblies, forming a single chromosomal contig. Nevertheless, assembly quality is influenced by read quality, coverage, and parameter selection, making the methodological choices essential.
In addition, the comparison to a reference genome after the assembly allows for the identification of nucleotide variants. Alignment tools optimized for long-read and reference comparisons, such as Minimap2, provide accurate genome-wide alignments [3]. Reference-based comparisons also facilitate validation of assembly accuracy and biological interpretation of observed differences. However, reliance on a reference genome may introduce bias if the sequenced strain differs substantially from the reference sequence.
Methods
Sequencing Data
Raw sequencing reads will be provided in FASTQ format, from Oxford Nanopore Technologies with R10 chemistry. The dataset has a read accuracy of Q20+ and a read length N50 between 5 and 15 kb, which is appropriate for long-read assembly of bacterial genomes.
Quality control
Initial quality assessment of raw reads will be performed to show and evaluatethe  read length distributions and overall dataset characteristics. Although long-read assemblers such as Flye can tolerate sequencing errors, maintaining sufficient read quality is important to reduce mis-assemblies and improve consensus accuracy [1]. Therefore, extremely low-quality reads should be excluded prior to assembly.
Genome assembly
De novo genome assembly will be conducted using Flye (v2.9 macos), a long-read assembler designed for error-prone sequencing data [2]. Flye will be run in --nano-raw mode, which is appropriate for uncorrected Oxford Nanopore reads. The genome size parameter will be set to approximately 5 Mb, reflecting the expected size of Salmonella enterica. This genome size estimate helps guide graph construction and repeat resolution during assembly [2].
Assembly post-processing and polishing
Post-assembly polishing will be performed to improve the consensus sequence by correcting potential system sequencing errors. Post-processing is particularly important for Nanopore-only assemblies, because residual base-calling errors may remain after initial assembly [1]. Polishing tools trained on Nanopore data will be used to refine the final consensus sequence.
Reference genome retrieval
A complete reference genome for Salmonella enterica will be downloaded from the NCBI RefSeq database, which provides standardized reference sequences with associated metadata [5]. The reference genome accession number and strain information will be documented to ensure reproducibility.
Genome alignment and comparison
The polished assembly will be aligned to the reference genome using Minimap2 (v2.26), which is optimized for fast and accurate alignment of long sequences and genome assemblies [3]. Default parameters suitable for assembly-to-reference alignment will be used. This alignment will serve as the basis for identifying sequence differences between the assembled genome and the reference.
Assembly quality check
Assembly quality metrics, including total assembly length, number of contigs, and N50, will be assessed using QUAST (v5.2.0) [4]. Comparing these metrics to the expected genome size of Salmonella enterica will help evaluate assembly completeness and structural accuracy.

References 
1.	Amarasinghe SL, Su S, Dong X, Zappia L, Ritchie ME, Gouil Q. Opportunities and challenges in long-read sequencing data analysis. Genome Biology. 2020;21:30.
2.	Kolmogorov M, Yuan J, Lin Y, Pevzner PA. Assembly of long, error-prone reads using repeat graphs. Nature Biotechnology. 2019;37(5):540–546.
3.	Li H. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics. 2018;34(18):3094–3100.
4.	Gurevich A, Saveliev V, Vyahhi N, Tesler G. QUAST: quality assessment tool for genome assemblies. Bioinformatics. 2013;29(8):1072–1075.
5.	O’Leary NA, Wright MW, Brister JR, et al. Reference sequence (RefSeq) database at NCBI: current status. Nucleic Acids Research. 2016;44(D1):D733–D745.

