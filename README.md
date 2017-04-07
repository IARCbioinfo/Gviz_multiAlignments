# Gviz_multiAlignments
R script for multiple BAM alignments viewing using Gviz (bioconductor package)

## Command line example:
```
Rscript script_gviz.r --pos_file=file_name.txt --bam_folder=/path_to_BAMs/ --ref_genome=Hsapiens.UCSC.hg19 --sample_names=SAMPLE
```
Arguments :

| Parameter | Default value | Description |
|-----------|--------------:|-------------|
| bam_folder    |            - | Folder containing all BAMs needed |
| pos_file | - |  Input file |
| ref_genome | - | String defining the reference genome used for alignment |
| sample_names | FILE | Set this argument to "SAMPLE" if the input file contains the samples names extracted from the BAM files and not the BAM files names |

#### Example of an input file (pos_file argument)
| | | | |
|-----------|--------------|-------------|-----------|
| chr17 | 7572814	| SampleName1 or BAM_file1 |	 |
| chr17	| 7572814	| SampleName1 or BAM_file1	| SampleName2	or BAM_file2 |
| chr17	| 7579643	| SampleName1 or BAM_file1	| SampleName2	or BAM_file2 |

Using this input file, the R script will generate 3 pdf files, each pdf file containing the alignment of the BAM file(s) at each position.  

#### Bioconductor packages to install for plotting the alignments :

- The [Gviz](https://bioconductor.org/packages/release/bioc/html/Gviz.html) package
- A [BSgenome data package](https://bioconductor.org/packages/release/BiocViews.html#___BSgenome) to provide a full genome sequence. This sequence can be, for Homo sapiens, provided by UCSC or based on NCBI GRCh37 for the 1000genomes Reference Genome Sequence ([hs37d5](https://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.1000genomes.hs37d5.html)). 
- An [Annotation package for TxDb objects](http://bioconductor.org/packages/release/BiocViews.html#___TxDb).
- A [Genome wide annotation](https://bioconductor.org/packages/release/BiocViews.html#___OrgDb), it contains mappings between Entrez Gene identifiers and GenBank accession numbers. Examples : the Genome wide annotation package for Human : *org.Hs.eg.db* and for the mouse : *org.Mm.eg.db*.

Examples :

If one needs the UCSC version of the reference Human genome hg19, the fowolling packages should be installed : 
- *Gviz*
- *BSgenome.Hsapiens.UCSC.hg19* (hg18 and hg38 UCSC version can also be used)
- *TxDb.Hsapiens.UCSC.hg19.knownGene* 
- *org.Hs.eg.db* 

These packages exist for other organisms than Human but have not been tested.

One can for example generate the alignments plot for data issued from the mouse by installing : 
- *BSgenome.Mmusculus.UCSC.mm10* (mm9 UCSC version can also be used)
- *TxDb.Mmusculus.UCSC.mm10.knownGene*  
- *org.Mm.eg.db*

For the other organisms the packages needs to have the same nomenclature as the ones listed above.

#### --ref_genome argument :

The name of the reference genome as to be the same as the BSgenome data package name without "BSgenome.", for the UCSC version of the reference Human genome hg19, one needs to set --ref_genome to "*Hsapiens.UCSC.hg19*".


#### The alignment plot from top to bottom :

- Chromosome representation : a red vertical line shows the position of the variant
- Genomic axis associated with the alignment.
- BAM alignment(s) (50 bases on both sides of the variant) : the variant position is highlighted in red.
- The reference genome : the variant position is highlighted in red.
- Genome annotation : the yellow blocks represent exons, the variant position is highlighted in red. The annotation is represented only if the variant is not in an intergenic region.
- Zoom-out on the genome annotation : representation of the whole gene whose name is on the right side of the gene annotation. The annotation is represented only if the variant is not in an intergenic region. A red vertical line shows the position of the variant. 
- Genomic axis corresponding to the previous genome annotation. A red vertical line shows the position of the variant. 
