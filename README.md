##Build instructions

###Packages Required

* gcc >= 4.4
* cmake
* zlib
* boost libraries >= 1.44
* google-sparse-hash

Note: it is advised to have installed the latest version of the previous packages.

###How to build the application

    $ cd gam-ngs
    $ mkdir build
    $ cd build
    $ cmake ..
    $ make

GAM-NGS's executables are put in the "bin" sub-folder.

##Usage

###Prerequisites

GAM-NGS needs in input, for each assembly and for each read library, a file that lists BAM files of aligned libraries.
This file must be formatted as shown in the following example:

    $ cat assembly.PE.bams.txt
    /path/to/bam/file/master-pe-lib1.bam
    <min-insert-size> <max-insert-size>
    /path/to/bam/file/master-pe-lib2.bam
    <min-insert-size> <max-insert-size>

where each bam's path is followed by a line that specifies the minimum and maximum insert size that should be used by GAM-NGS to compute insert size mean and standard deviation.
Moreover, all provided BAM files have to be coordinate-sorted (command samtools sort \<in.bam\> \<out.prefix\>), along with the corresponding index file (command: samtools index \<in.sorted.bam\>).

###Blocks' construction

    $ gam-create --master-bam <master.PE.bams.txt> --slave-bam <slave.PE.bams.txt> --min-block-size <min-reads> --output <output.prefix>

where \<min-reads\> is the number of reads required to build a block (region with the same reads aligned in master/slave assemblies).

The previous command will create the following files:
- \<output.prefix\>.blocks        blocks descriptor
- \<master.PE.bams.txt\>.isize    libraries' statistics (insert size mean, standard deviation, read coverage)
- \<slave.PE.bams.txt\>.isize     libraries' statistics (insert size mean, standard deviation, read coverage)

###Merging

    $gam-merge --master-bam <master.PE.bams.txt> --slave-bam <slave.PE.bams.txt> --blocks-file <blocks-file> --master-fasta <master.fasta> --slave-fasta <slave.fasta> --min-block-size <min-block-size> --threads <threads> --output <output.prefix> 2> merge.err

where:
* \<blocks-file\> is the blocks descriptor file created with gam-create command.
* \<min-block-size\> specifies the minimum number of reads a block must have to be used.
* \<threads\> specifies the number of threads used in the merging phase.

The previous command will create the following files:
* \<output.prefix\>.gam.fasta            merged assembly
* \<output.prefix\>.pctgs                merged contigs descriptor
* \<output.prefix\>.noblocks.BF.fasta    slave contigs without blocks (before GAM-NGS's filtering phase)
* \<output.prefix\>.noblocks.AF.fasta    slave contigs without blocks (after GAM-NGS's filtering phase)
* \<master.PE.bams.txt\>.isize           libraries' statistics (if not previously created with gam-create)
* \<slave.PE.bams.txt\>.isize            libraries' statistics (if not previously created with gam-create)


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


##Bug reporting

If gdb package is available in your system and you found a bug in GAM-NGS (e.g., segmentation fault),
I will appreciate if you could run the tool in "debug mode" and send me the log file.
This would help me a lot in finding the source of the problem.

The steps are the following ones:

 1. Edit CMakeLists.txt file at line 17, replacing "Release" with "Debug" so that the line is:

	set(CMAKE_BUILD_TYPE Debug)

 2. Re-compile the tool with the usual commands (removing any previous file/folder in the build directory)
	
	$ cmake ..
	$ make

 3. Run gdb with the executable (gam-create or gam-merge) in which the problem is present:

	$ gdb /path/to/gam-create-or-merge

 4. Run the executable specifying the arguments that caused the problem (possibly using a single thread):

	(gdb) set logging on
	(gdb) run --master-bam <master-pe-list.txt> --slave-bam <slave-pe-list.txt> ...

 5. After the problem has occurred, please run the following command and exit gdb:

	(gdb) backtrace
	(gdb) quit

 6. A file named gdb.txt should have been created in the current directory. 
 Attach that file to help me solve the problem (possibly along with standard output and standard error).