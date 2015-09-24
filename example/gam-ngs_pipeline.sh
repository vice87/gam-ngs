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
BLOCKSIZE=10
BWA=bwa

master=Allpaths-LG
slave=MSR-CA

GAM_CREATE=../bin/gam-create
GAM_MERGE=../bin/gam-merge

# Aligning reads on master and slave's assemblies

>./log.stderr
for index in $master $slave
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

gam-ngs () {
	# Variable assignment
	MASTER=$1
	SLAVE=$2
	# GAM-NGS pipeline

	mkdir -p ${SCRIPT_PATH}/gam-ngs_merge_${MASTER}_${SLAVE}

	# Prepare GAM-NGS's input files
	for name in $MASTER $SLAVE
	do
	        echo -e "${SCRIPT_PATH}/Alignments/${name}/${name}.frag.bwa.sorted.bam\n90 270" >${SCRIPT_PATH}/gam-ngs_merge_${MASTER}_${SLAVE}/${name}.PE.list.txt
	done

	# GAM-NGS: Blocks construction phase

	echo -e "\n### GAM-NGS's blocks construction\n"
	command="${GAM_CREATE} --master-bam ${SCRIPT_PATH}/gam-ngs_merge_${MASTER}_${SLAVE}/${MASTER}.PE.list.txt --slave-bam ${SCRIPT_PATH}/gam-ngs_merge_${MASTER}_${SLAVE}/${SLAVE}.PE.list.txt --min-block-size ${BLOCKSIZE} --output ${SCRIPT_PATH}/gam-ngs_merge_${MASTER}_${SLAVE}/out >${SCRIPT_PATH}/gam-ngs_merge_${MASTER}_${SLAVE}/gam-create.log.out 2>${SCRIPT_PATH}/gam-ngs_merge_${MASTER}_${SLAVE}/gam-create.log.err"

	echo $command
	eval $command

	# GAM-NGS: Merge

	echo -e "### GAM-NGS's merging phase\n"
	command="${GAM_MERGE} --blocks-file ${SCRIPT_PATH}/gam-ngs_merge_${MASTER}_${SLAVE}/out.blocks --master-bam ${SCRIPT_PATH}/gam-ngs_merge_${MASTER}_${SLAVE}/${MASTER}.PE.list.txt --master-fasta ${SCRIPT_PATH}/Assembly/${MASTER}/genome.ctg.fasta --slave-bam ${SCRIPT_PATH}/gam-ngs_merge_${MASTER}_${SLAVE}/${SLAVE}.PE.list.txt --slave-fasta ${SCRIPT_PATH}/Assembly/${SLAVE}/genome.ctg.fasta --min-block-size ${BLOCKSIZE} --output ${SCRIPT_PATH}/gam-ngs_merge_${MASTER}_${SLAVE}/out --threads ${THREADS_NUM} >${SCRIPT_PATH}/gam-ngs_merge_${MASTER}_${SLAVE}/gam-merge.log.out 2>${SCRIPT_PATH}/gam-ngs_merge_${MASTER}_${SLAVE}/gam-merge.log.err"

	echo $command
	eval $command
}
# First run
gam-ngs $master $slave

# Second run
gam-ngs $slave $master


