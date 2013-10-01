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
