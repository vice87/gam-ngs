#!/bin/bash

#  This file is part of GAM-NGS.
#  Copyright (c) 2011 by Riccardo Vicedomini <rvicedomini@appliedgenomics.org>,
#  Francesco Vezzi <vezzi@appliedgenomics.org>,
#  Simone Scalabrin <scalabrin@appliedgenomics.org>,
#  Lars Arverstad <lars.arvestad@scilifelab.se>,
#  Alberto Policriti <policriti@appliedgenomics.org>,
#  Alberto Casagrande <casagrande@appliedgenomics.org>
#
#  GAM-NGS is an evolution of a previous work (GAM) done by Alberto Casagrande,
#  Cristian Del Fabbro, Simone Scalabrin, and Alberto Policriti.
#  In particular, GAM-NGS has been adapted to work on NGS data sets and it has
#  been written using GAM's software as starting point. Thus, it shares part of
#  GAM's source code.
#
#  Moreover, GAM-NGS uses BamTools library to access BAM files.
#  BamTools's source code has been put in ./lib/bamtools-2.0.5/ folder,
#  in which its license can be found.
#
#  GAM-NGS is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  GAM-NGS is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with GAM-NGS.  If not, see <http://www.gnu.org/licenses/>.

SCRIPT_PATH=`pwd -P`
THREADS_NUM=4
BWA=bwa
GAM_CREATE=../bin/gam-create
GAM_MERGE=../bin/gam-merge

# Aligning reads on Allpaths-LG and MSR-CA's assemblies

>./log.stderr

for index in Allpaths-LG MSR-CA
do
	ASM_PATH="./Assembly/${index}/genome.ctg.fasta"
	ALN_PREFIX="./Alignments/${index}/${index}.frag.bwa"
	PE1_PATH="./Data/allpathsCor/frag_1.fastq"
	PE2_PATH="./Data/allpathsCor/frag_2.fastq"
	SAI1_PATH="./Alignments/${index}/frag_1.sai"
	SAI2_PATH="./Alignments/${index}/frag_2.sai"

	echo -e "\n### Aligning PE reads on ${index} assembly"

	echo -e "\nmkdir -p ./Alignments/${index}"
	mkdir -p ./Alignments/${index}

	echo "bwa index $ASM_PATH"
	$BWA index $ASM_PATH 2>>./log.stderr

	echo "bwa aln -o 0 -t $THREADS_NUM $ASM_PATH $PE1_PATH > $SAI1_PATH"
	$BWA aln -o 0 -t $THREADS_NUM $ASM_PATH $PE1_PATH > $SAI1_PATH 2>>./log.stderr

	echo "bwa aln -o 0 -t $THREADS_NUM $ASM_PATH $PE2_PATH > $SAI2_PATH"
	$BWA aln -o 0 -t $THREADS_NUM $ASM_PATH $PE2_PATH > $SAI2_PATH 2>>./log.stderr

	echo "bwa sampe $ASM_PATH $SAI1_PATH $SAI2_PATH $PE1_PATH $PE2_PATH | samtools view -Shu -@ $THREADS_NUM - | samtools sort -@ $THREADS_NUM - $ALN_PREFIX.sorted"
	($BWA sampe $ASM_PATH $SAI1_PATH $SAI2_PATH $PE1_PATH $PE2_PATH | samtools view -Sb -@ $THREADS_NUM - | samtools sort -@ $THREADS_NUM - $ALN_PREFIX.sorted) 2>>./log.stderr

	rm $SAI1_PATH
	rm $SAI2_PATH

	echo "samtools index $ALN_PREFIX.sorted.bam"
	samtools index $ALN_PREFIX.sorted.bam 2>>./log.stderr
done

# GAM-NGS pipeline

mkdir -p ${SCRIPT_PATH}/gam-ngs_merge

# Prepare GAM-NGS's input files

echo -e "${SCRIPT_PATH}/Alignments/Allpaths-LG/Allpaths-LG.frag.bwa.sorted.bam\n90 270" >${SCRIPT_PATH}/gam-ngs_merge/Allpaths-LG.PE.list.txt
echo -e "${SCRIPT_PATH}/Alignments/MSR-CA/MSR-CA.frag.bwa.sorted.bam\n90 270" >${SCRIPT_PATH}/gam-ngs_merge/MSR-CA.PE.list.txt

# GAM-NGS: Blocks construction phase

echo -e "\n### GAM-NGS's blocks construction\n"
echo -e "gam-create --master-bam ${SCRIPT_PATH}/gam-ngs_merge/Allpaths-LG.PE.list.txt --slave-bam ${SCRIPT_PATH}/gam-ngs_merge/MSR-CA.PE.list.txt --min-block-size 10 --output ${SCRIPT_PATH}/gam-ngs_merge/out >${SCRIPT_PATH}/gam-ngs_merge/gam-create.log.out 2>${SCRIPT_PATH}/gam-ngs_merge/gam-create.log.err\n"
${GAM_CREATE} --master-bam ${SCRIPT_PATH}/gam-ngs_merge/Allpaths-LG.PE.list.txt --slave-bam ${SCRIPT_PATH}/gam-ngs_merge/MSR-CA.PE.list.txt --min-block-size 10 --output ${SCRIPT_PATH}/gam-ngs_merge/out >${SCRIPT_PATH}/gam-ngs_merge/gam-create.log.out 2>${SCRIPT_PATH}/gam-ngs_merge/gam-create.log.err

# GAM-NGS: Merge

echo -e "### GAM-NGS's merging phase\n"
echo -e "gam-merge --blocks-file ${SCRIPT_PATH}/gam-ngs_merge/out.blocks --master-bam ${SCRIPT_PATH}/gam-ngs_merge/Allpaths-LG.PE.list.txt --master-fasta ${SCRIPT_PATH}/Assembly/Allpaths-LG/genome.ctg.fasta --slave-bam ${SCRIPT_PATH}/gam-ngs_merge/MSR-CA.PE.list.txt --slave-fasta ${SCRIPT_PATH}/Assembly/MSR-CA/genome.ctg.fasta --min-block-size 10 --output ${SCRIPT_PATH}/gam-ngs_merge/out --threads ${THREADS_NUM} >${SCRIPT_PATH}/gam-ngs_merge/gam-merge.log.out 2>${SCRIPT_PATH}/gam-ngs_merge/gam-merge.log.err\n"
${GAM_MERGE} --blocks-file ${SCRIPT_PATH}/gam-ngs_merge/out.blocks --master-bam ${SCRIPT_PATH}/gam-ngs_merge/Allpaths-LG.PE.list.txt --master-fasta ${SCRIPT_PATH}/Assembly/Allpaths-LG/genome.ctg.fasta --slave-bam ${SCRIPT_PATH}/gam-ngs_merge/MSR-CA.PE.list.txt --slave-fasta ${SCRIPT_PATH}/Assembly/MSR-CA/genome.ctg.fasta --min-block-size 10 --output ${SCRIPT_PATH}/gam-ngs_merge/out --threads ${THREADS_NUM} >${SCRIPT_PATH}/gam-ngs_merge/gam-merge.log.out 2>${SCRIPT_PATH}/gam-ngs_merge/gam-merge.log.err
