# Continuous Analysis RNA-Seq Differential Expression Analysis Example

This is a sample repository showing a [Continuous Analysis Workflow](https://github.com/greenelab/continuous_analysis) for RNA-Seq analysis. A description of continuous analysis is available as a [pre-print](http://biorxiv.org/content/early/2016/06/01/056473).

In this example we follow the workflow described by [David Balli](https://benchtobioinformatics.wordpress.com/2015/07/10/using-kallisto-for-gene-expression-analysis-of-published-rnaseq-data/) and use data generated from [Boj et al.](http://www.cell.com/cell/abstract/S0092-8674(14)01592-X) (open access). Balli used a similar workflow to the one described by [Andrew Mckenzie](https://andrewtmckenzie.com/2015/05/12/how-to-run-kallisto-on-ncbi-sra-rna-seq-data-for-differential-expression-using-the-mac-terminal/).

To preform this analysis we use several tools:

* [Kallisto](https://pachterlab.github.io/kallisto/) 
* [Limma](https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf)
* [Sleuth](http://pachterlab.github.io/sleuth/)

A variant of this Continuous Analysis Workflow that performs the same analysis, but which uses [Salmon](https://combine-lab.github.io/salmon) for transcript quantification is available [here](https://github.com/COMBINE-lab/continuous_analysis_rnaseq).  This analysis also demonstrates the use of [Drone v0.5](https://github.com/drone/drone), which eliminates the need to fetch the data into the image on each run.

## Sample Results

The Continuous Analysis process generates several useful artifacts, two of which are (the full experiment is described below):

1. Change logs/synchronization between code and figures.

![](https://raw.githubusercontent.com/greenelab/continuous_analysis_rnaseq/master/references/pca.png)
The effect of adding an additional organoid derived from pancreatic adenocarcinoma on principal components analysis using Kallisto’s estimated counts.

![](https://raw.githubusercontent.com/greenelab/continuous_analysis_rnaseq/master/references/volcano.png)
A volcano plot plotting the p-value vs. the log fold change. Adding an additional organoid derived from pancreatic adenocarcinoma leads to an additional gene being marked as significantly differentially expressed after Benjamini & Hochberg correction.


2. Complete "audit" logs of the code run: [Logs](https://raw.githubusercontent.com/greenelab/continuous_analysis_rnaseq/master/references/full_logs.txt)

![](https://raw.githubusercontent.com/greenelab/continuous_analysis_rnaseq/master/references/logs.png)



## Description of analysis
We followed a reduced analysis workflow demonstrated by Balli using the SRA files for 8 samples: 2 normal, 3 mP, 3 mT). These samples represent extract to approximately 480 million reads and 150gb of data (FASTQ format). We perform this experiment first with 7 samples (2 normal, 3mP, 2mT) and then add the 8th sample to view the differences.

We perform two preprocessing steps prior to beginning continuous analysis (details/reasoning in continuous analysis configuration section).

Download the samples from the [Sequence Read Archive](http://www.ncbi.nlm.nih.gov/sra?term=SRP049959 "SRR1654626", "SRR1654628", "SRR1654633", "SRR1654636", "SRR16546367", “SRR1654639”, "SRR1654637", "SRR1654641", "SRR1654643")
Split the .sra into fast q files using the [SRA toolkit](http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc)
Download the [mouse reference genome assembly](http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/refMrna.fa.gz) 

Continuous Analysis Run ([script](https://github.com/greenelab/continuous_analysis_rnaseq/blob/master/.drone.yml)):

1. Generate a kallisto index file from the reference file and quantify abundances of transcripts from each RNA-Seq sample (run on 28 cores).
2. The next portion of the analysis is performed from ‘r_script.r’ and follows the workflow described by Balli: Generate the transcripts per million (TPM) matrix.
3. Create a matrix to specify the group each sample belongs to.
4. Filter out lowly expressed genes.
5. Generate a principle component plot
6. Fit the limma linear model for differential gene expression analysis.
7. Plot differential expression in the form of a volcano plot.


## Feedback

Please feel free to email me - (brettbe) at med.upenn.edu with any feedback or raise a github issue with any comments or questions.

## Acknowledgements

We would like to thank David Balli for his post providing the analysis design and significant source code used in this example.

This work is supported by the Gordon and Betty Moore Foundation's Data-Driven Discovery Initiative through Grant GBMF4552 to C.S.G. as well as the Commonwealth Universal Research Enhancement (CURE) Program grant from the Pennsylvania Department of Health.
