#!/bin/bash
module load bwa
module load openjdk
module load bioawk
module load pigz

# setting up variables
read1=$(realpath $1)
read2=$(realpath $2)
refSeq=$(realpath $3)
projectDir="/ptmp/arnstrm/hic_mapping/hap1_test_run"

threads=${SLURM_CPUS_ON_NODE}
splitthreads=$(expr ${threads} / 2)
memory=$(free -g | awk 'NR==2{print $4"G"}')

juicer_version=1.6.2
ext=$(basename ${projectDir})
site="DpnII"
ligation="GATCGATC"
groupname="a$(date +%s)"
justexact=0
genomeID="${ext}"
splitsize=900000000000
baseoutname=$(basename ${read1} |cut -f1 -d "_")


# get scripts for the pipeline
cd ${projectDir}
git clone git@github.com:HuffordLab/Juicer_PanAnd.git
mv Juicer_PanAnd juicerdir
juiceDir="${projectDir}/juicerdir"

# directory structure:
mkdir -p ${projectDir}/fastq
mkdir -p ${projectDir}/references
mkdir -p ${projectDir}/aligned
mkdir -p ${projectDir}/splits
mkdir -p ${projectDir}/debug
mkdir -p ${projectDir}/restriction_sites
mkdir -p ${projectDir}/HIC_tmp

outputdir="${projectDir}/aligned"
splitdir="${projectDir}/splits"
debugdir="${projectDir}/debug"
tmpdir="${projectDir}/HIC_tmp"

# files in right places
ln -s ${read1} ${projectDir}/fastq/
ln -s ${read2} ${projectDir}/fastq/
ln -s ${refSeq} ${projectDir}/references/

# not sure if this is necessary, but it's in the original script
# pigz -p ${threads} -d -c ${read1} > ${splitdir}/${baseoutname}_R1.fastq000.fastq
# pigz -p ${threads} -d -c ${read2} > ${splitdir}/${baseoutname}_R2.fastq000.fastq
# probably can get away with just dummy files in there instead
# touch ${splitdir}/${baseoutname}_R1.fastq000.fastq
# touch ${splitdir}/${baseoutname}_R2.fastq000.fastq

# chrom sizes file
bioawk -c fastx '{print $name"\t"length($seq)}' ${refSeq} > ${projectDir}/chrom.sizes
genomePath=${projectDir}/chrom.sizes

# restriction sites for the genome
python ${juiceDir}/scripts/generate_site_positions.py DpnII ${refSeq%.*} ${refSeq}
mv ${refSeq%.*}_${site}.txt ${projectDir}/restriction_sites/
site_file="${projectDir}/restriction_sites/${refSeq%.*}_${site}.txt"

# genome index
cd ${projectDir}/references
bwa index ${refSeq} && cd ${projectDir}
refSeq="${projectDir}/references/$(basename ${refSeq})"

# creating header file
headfile="${projectDir}/debug/headers"
echo -e "Juicer version $juicer_version;"  > ${headfile}
echo -e "BWA: $(bwa 2>&1 | awk '$1 ~/Version/ {print $NF}')" >> ${headfile}
echo -e "$threads threads;\n${memory} memory" >> ${headfile}
java --version | head -n 1 >> ${headfile}
${juiceDir}/scripts/juicer_tools -V | head -n 1 >> ${headfile}
echo -e "juicer.sh -d ${projectDir} -s ${site} -p ${genomePath} -y ${site_file} -z ${refSeq} -C ${splitsize} -D ${juicerDir} -b ${ligation} -t ${threads}" >> ${headfile}

# countLigations:
num1=$(paste <(bioawk -c fastx '{print $seq}' ${read1}) <(bioawk -c fastx '{print $seq}' ${read2}) | grep -cE ${ligation})
num2=$(pigz -d -p ${threads} -c ${read1} | wc -l | awk '{print $1}')
echo -e "${num1}" > ${splitdir}/${ext}_${baseoutname}_norm.txt.res.txt
echo -e "${num2}" > ${splitdir}/${ext}_${baseoutname}_linecount.txt

# bwa align
bwa mem -SP5M -t $SLURM_JOB_CPUS_PER_NODE ${refSeq} ${read1} ${read2} > ${splitdir}/${ext}_${baseoutname}.sam

# chimera_blacklist
awk -v "fname1"=${splitdir}/${ext}_${baseoutname}_norm.txt -v "fname2"=${splitdir}/${ext}_${baseoutname}_abnorm.sam -v "fname3"=${splitdir}/${ext}_${baseoutname}_unmapped.sam -f ${juiceDir}/scripts/chimeric_blacklist.awk ${splitdir}/${ext}_${baseoutname}.sam

# fragment assignment
perl ${juiceDir}/scripts/fragment.pl ${splitdir}/${ext}_${baseoutname}_norm.txt ${splitdir}/${ext}_${baseoutname}.frag.txt ${site_file}

# sort by chr, frag, strand, and position
sort --parallel=${threads} \
     --buffer-size=${memory} \
     --temporary-directory=${tmpdir} \
     --key=2,2d --key=6,6d --key=4,4n --key=8,8n --key=1,1n --key=5,5n --key=3,3n ${splitdir}/${ext}_${baseoutname}.frag.txt > ${splitdir}/${ext}_${baseoutname}.sort.txt

# fragmerge
cp ${splitdir}/${ext}_${baseoutname}.sort.txt ${outputdir}/merged_sort.txt

# split_rmdups
awk -v groupname=${groupname} -v debugdir=${debugdir} -v juicedir=${juiceDir} -v site=${site} -v genomeID=${genomeID} -v genomePath=${genomePath} -v justexact=0 -f ${juiceDir}/scripts/split_rmdups.awk ${outputdir}/merged_sort.txt

# prestats
tail -n1 ${headfile} | awk '{printf"%-1000s\n", $0}' > ${outputdir}/inter.txt
cat ${splitdir}/*.res.txt |\
    awk -f ${juiceDir}/scripts/stats_sub.awk >> ${outputdir}/inter.txt
${juiceDir}/scripts/juicer_tools LibraryComplexity ${outputdir} inter.txt >> ${outputdir}/inter.txt
cp ${outputdir}/inter.txt ${outputdir}/inter_30.txt


perl ${juiceDir}/scripts/statistics.pl \
    -s ${site_file} \
    -l ${ligation} \
    -o ${outputdir}/inter.txt \
    -q 1 ${outputdir}/merged_nodups.txt

perl ${juiceDir}/scripts/statistics.pl \
    -s ${site_file} \
    -l ${ligation} \
    -o ${outputdir}/inter_30.txt \
    -q 30 ${outputdir}/merged_nodups.txt

#collect collisions and dedup
cat $splitdir/*_abnorm.sam > $outputdir/abnormal.sam
cat $splitdir/*_unmapped.sam > $outputdir/unmapped.sam
awk -f ${juiceDir}/scripts/collisions.awk ${outputdir}/abnormal.sam > ${outputdir}/collisions.txt
gawk -v fname=${outputdir}/collisions.txt -f ${juiceDir}/scripts/collisions_dedup_rearrange_cols.awk ${outputdir}/collisions.txt |\
      sort -k3,3n -k4,4n -k10,10n -k11,11n -k17,17n -k18,18n -k24,24n -k25,25n -k31,31n -k32,32n |\
      awk -v name=${outputdir} -f ${juiceDir}/scripts/collisions_dups.awk

#${juiceDir}/scripts/juicer_tools pre -s /${outputdir}/inter.txt -g ${outputdir}/inter_hists.m -q 1 ${outputdir}/merged_nodups.txt ${outputdir}/inter.hic ${genomePath}
${juiceDir}/scripts/juicer_tools pre -f ${site_file} -s ${outputdir}/inter.txt -g ${outputdir}/inter_hists.m -q 1 ${outputdir}/merged_nodups.txt ${outputdir}/inter.hic ${genomePath}

#${juiceDir}/scripts/juicer_tools pre -s ${outputdir}/inter_30.txt -g ${outputdir}/inter_30_hists.m -q 30 ${outputdir}/merged_nodups.txt ${outputdir}/inter_30.hic ${genomePath}
${juiceDir}/scripts/juicer_tools pre -f ${site_file} -s ${outputdir}/inter_30.txt -g ${outputdir}/inter_30_hists.m -q 30 ${outputdir}/merged_nodups.txt ${outputdir}/inter_30.hic ${genomePath}

${juiceDir}/scripts/juicer_arrowhead.sh -j ${juiceDir}/scripts/juicer_tools -i ${outputdir}/inter_30.hic

${juiceDir}/scripts/check.sh
