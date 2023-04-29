#!/bin/bash

# setting up variables
read1=$1
read2=$2
refSeq=$3
ext="tdats"
juiceDir="/ptmp/arnstrm/hic_mapping/hap1_test_run/juicerdir"
topDir=${juiceDir}
juicer_version=1.6.2

threads=${SLURM_CPUS_ON_NODE}
memory=$(free -g | awk 'NR==2{print $4"G"}')
name=$(basename ${read1%%.*})
site="DpnII"
ligation="GATCGATC"
queue="amd"
long_queue=${queue}
headfile="${juiceDir}/debug/headers"
groupname="a$(date +%s)"
justexact=0

genomeID="${ext}"
genomePath=${juiceDir}/chrom.sizes

splitsize=900000000000

mkdir -p ${juiceDir}/references
mkdir -p ${juiceDir}/aligned
mkdir -p ${juiceDir}/splits
mkdir -p ${juiceDir}/debug
mkdir -p ${juiceDir}/restriction_sites
mkdir -p ${juiceDir}/HIC_tmp
# restriction sites for the genome
python ${juiceDir}/scripts/generate_site_positions.py DpnII ${refSeq%.*} ${refSeq}
mv ${refSeq%.*}_${site}.txt ${juiceDir}/restriction_sites/
site_file="${juiceDir}/restriction_sites/${refSeq%.*}_${site}.txt"


outputdir="${juiceDir}/aligned"
splitdir="${juiceDir}/splits"
debugdir="${juiceDir}/debug"
tmpdir="${juiceDir}/HIC_tmp"

#need to change this to be more flexible
baseoutname=$(basename ${read1} |cut -f1 -d "_")
guardjid=""
sortthreadstring="-t ${threads}

# header script
cat > head.sh << EOF
module load bwa
module load openjdk
echo -e "Juicer version $juicer_version;" 
echo -e "BWA: $(bwa 2>&1 | awk '$1 ~/Version/ {print $NF}')"
echo -e "$threads threads; "
java --version |head -n 1
${juiceDir}/scripts/juicer_tools -V | head -n 1
echo "$0 $@"
EOF
chmod + x head.sh
head.sh > ${headfile}

# countLigations:
num1=$(paste <(gunzip -c ${read1}) <(gunzip -c ${read2}) | awk '!((NR+2)%4)' | grep -cE ${ligation})
num2=$(gunzip -c ${read1} | wc -l | awk '{print $1}')
echo -e "${num1}" > ${splitdir}/${ext}_${baseoutname}_norm.txt.res.txt
echo -e "${num2}" > ${splitdir}/${ext}_${baseoutname}_linecount.txt

# bwa align
bwa mem -SP5M -t $SLURM_JOB_CPUS_PER_NODE ${refSeq} ${read1} ${read2} > ${splitdir}/${ext}_${baseoutname}.sam

# chimera_blacklist
awk -v "fname1"=${splitdir}/${ext}_${baseoutname}_norm.txt \
    -v "fname2"=${splitdir}/${ext}_${baseoutname}_abnorm.sam \
    -v "fname3"=${splitdir}/${ext}_${baseoutname}_unmapped.sam \
    -f ${juiceDir}/scripts/chimeric_blacklist.awk ${splitdir}/${ext}_${baseoutname}.sam

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
awk -v queue="amd" \ 
    -v groupname=${groupname} \
    -v debugdir=${debugdir} \
    -v dir=${outputdir} \
    -v topDir=${topDir} \
    -v juicedir=${juiceDir} \
    -v site=${site} \
    -v genomeID=${genomeID} \
    -v genomePath=${genomePath} \
    -v user=las \
    -v guardjid=4450589 \
    -v justexact=0 -f ${juiceDir}/scripts/split_rmdups.awk ${outputdir}/merged_sort.txt


# prestats
tail -n1 ${headfile} | awk '{printf"%-1000s\n", \$0}' > ${outputdir}/inter.txt
cat ${splitdir}/*.res.txt | awk -f ${juiceDir}/scripts/stats_sub.awk >> ${outputdir}/inter.txt
${juiceDir}/scripts/juicer_tools LibraryComplexity ${outputdir} inter.txt >> ${outputdir}/inter.txt
cp ${outputdir}/inter.txt ${outputdir}/inter_30.txt


perl ${juiceDir}/scripts/statistics.pl -s ${site_file} -l ${ligation} -o ${outputdir}/inter.txt -q 1 ${outputdir}/merged_nodups.txt
perl ${juiceDir}/scripts/statistics.pl -s ${site_file} -l ${ligation} -o ${outputdir}/inter_30.txt -q 30 ${outputdir}/merged_nodups.txt

#collect collisions and dedup
awk -f ${juiceDir}/scripts/collisions.awk ${outputdir}/abnormal.sam > ${outputdir}/collisions.txt
gawk -v fname=${outputdir}/collisions.txt \
     -f ${juiceDir}/scripts/collisions_dedup_rearrange_cols.awk ${outputdir}/collisions.txt |\
      sort -k3,3n -k4,4n -k10,10n -k11,11n -k17,17n -k18,18n -k24,24n -k25,25n -k31,31n -k32,32n |\
      awk -v name=${outputdir} -f ${juiceDir}/scripts/collisions_dups.awk

#${juiceDir}/scripts/juicer_tools pre -s /${outputdir}/inter.txt -g ${outputdir}/inter_hists.m -q 1 ${outputdir}/merged_nodups.txt ${outputdir}/inter.hic ${genomePath}
${juiceDir}/scripts/juicer_tools pre -f ${site_file} -s ${outputdir}/inter.txt -g ${outputdir}/inter_hists.m -q 1 ${outputdir}/merged_nodups.txt ${outputdir}/inter.hic ${genomePath}

#${juiceDir}/scripts/juicer_tools pre -s ${outputdir}/inter_30.txt -g ${outputdir}/inter_30_hists.m -q 30 ${outputdir}/merged_nodups.txt ${outputdir}/inter_30.hic ${genomePath}
${juiceDir}/scripts/juicer_tools pre -f ${site_file} -s ${outputdir}/inter_30.txt -g ${outputdir}/inter_30_hists.m -q 30 ${outputdir}/merged_nodups.txt ${outputdir}/inter_30.hic ${genomePath}

${juiceDir}/scripts/juicer_arrowhead.sh -j ${juiceDir}/scripts/juicer_tools -i ${outputdir}/inter_30.hic

${juiceDir}/scripts/check.sh
