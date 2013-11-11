##Example

### Prerequisites

The following packages are needed to run the example pipeline:

 * bwa: this example has been tested with version 0.7.5a-r405
 * samtools: this example has been tested with version 0.1.19-44428cd

### Download the data set

From the example directory execute the following commands:

    $ chmod +x ./download_dataset.sh
    $ ./download_dataset.sh

These previous commands will download Staphylococcus data set from GAGE. 
In particular, after the download, the following files will be used in GAM-NGS's usage example:

 * ./Data/allpathsCor/frag_1.fastq
 * ./Data/allpathsCor/frag_2.fastq
 * ./Assembly/Allpaths-LG/genome.ctg.fasta
 * ./Assembly/MSR-CA/genome.ctg.fasta

### Run the example pipeline

From the example directory execute the following commands (after downloading the data set):

    $ chmod +x ./gam-ngs_pipeline.sh
    $ ./gam-ngs_pipeline.sh

These will create alignment files (BAM) in ./Alignments sub-folder.
Then the merging with GAM-NGS of Allpaths-LG and MSR-CA assemblies will be performed in ./gam-ngs_merge sub-folder.