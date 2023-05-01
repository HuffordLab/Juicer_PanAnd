# Running Juicer manually

Request an interactive node, connected via `sleep` or `tmux` so that you don't lose progress when disconnected.

Following programs are required:

1. `bwa`
2. `samtools`
3. `bioawk`
4. `pigz`


## Step 1. Setting up project

#### 1. Input variables

```bash
read1="/ptmp/arnstrm/hic_mapping/AP8220001_R1.fastq.gz"
read2="/ptmp/arnstrm/hic_mapping/AP8220001_R2.fastq.gz"
refSeq="/ptmp/arnstrm/hic_mapping/tdacts_asm_hap1.fasta"
projectDir="/ptmp/arnstrm/hic_mapping/tdacts_asm_hap1"
```

#### 2. Runtime variables

```bash
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
```

#### 3. Scripts folder

This is main Juicer version 1.6.

```bash
cd ${projectDir}
git clone git@github.com:HuffordLab/Juicer_PanAnd.git
mv Juicer_PanAnd juicerdir
juiceDir="${projectDir}/juicerdir"
```

#### 4. Directory structure:

```bash
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
```

#### 5. Files in correct places

```bash
ln -s ${read1} ${projectDir}/fastq/
ln -s ${read2} ${projectDir}/fastq/
ln -s ${refSeq} ${projectDir}/references/
# touch ${splitdir}/${baseoutname}_R1.fastq000.fastq
# touch ${splitdir}/${baseoutname}_R2.fastq000.fastq
```


## Step 2. Preliminary files

Juicer requries genome index (BWA), restriction sites, and scaffold sizes  for the input genome. While the scaffold file size generation is fast, the other 2 files need considerable time depending on your genome size.

```bash
# chrom sizes file
bioawk -c fastx '{print $name"\t"length($seq)}' ${refSeq} > ${projectDir}/chrom.sizes
genomePath=${projectDir}/chrom.sizes
```

```bash
# restriction sites for the genome
python ${juiceDir}/scripts/generate_site_positions.py DpnII ${refSeq%.*} ${refSeq}
mv ${refSeq%.*}_${site}.txt ${projectDir}/restriction_sites/
site_file="${projectDir}/restriction_sites/$(basename ${refSeq%.*})_${site}.txt"
```

```bash
# bwa genome index
cd ${projectDir}/references
bwa index ${refSeq} && cd ${projectDir}
refSeq="${projectDir}/references/$(basename ${refSeq})"
```


## Step 3. Collect run metadata

Some information regarding the versions of programs is needed for the later summary stats files. Juicer captures the command you ran and saves it in the file along with the versions of various programs. We are going to simulate the same behaviour using the below commands. This is optional and is not needed for the generating the output file.


```bash
# creating header file
headfile="${projectDir}/debug/headers"
echo -en "Juicer version $juicer_version;"  > ${headfile}
echo -en "BWA: $(bwa 2>&1 | awk '$1 ~/Version/ {print $NF}')" >> ${headfile}
echo -en "$threads threads;\n${memory} memory" >> ${headfile}
java --version | head -n 1 >> ${headfile}
${juiceDir}/scripts/juicer_tools -V | head -n 1 >> ${headfile}
cat >> ${headfile} <<- EOF
juicer.sh \
-d ${projectDir} \
-s ${site} \
-p ${genomePath} \
-y ${site_file} \
-z ${refSeq} \
-D ${juicerDir} \
-b ${ligation} \
-t ${threads}
EOF
```

```bash
# countLigations
num1=$(paste <(bioawk -c fastx '{print $seq}' ${read1}) <(bioawk -c fastx '{print $seq}' ${read2}) | grep -cE ${ligation})
num2=$(pigz -d -p ${threads} -c ${read1} | wc -l | awk '{print $1}')
echo -en "${num1} " > ${splitdir}/${ext}_${baseoutname}_norm.txt.res.txt
echo -e "${num2}" > ${splitdir}/${ext}_${baseoutname}_linecount.txt
```

## Step 4. Processing

### bwa mapping

This probably takes a long time. Depening on how many CPUs you have and the memory. If possible, run this separately as a slurm job.

```bash
bwa mem -SP5M -t ${threads} ${refSeq} ${read1} ${read2} > ${splitdir}/${ext}_${baseoutname}.sam
```
For a 3Gb genome, with 40Gb x 2 compressed fastq files (517,952,838 x 2 reads), it took about 39,004.633 seconds walltime (10.8 hrs) or 2,463,509.733 seconds CPU time (684.31 hrs)

The machine:
```
Intel(R) Xeon(R) Gold 6140 CPU @ 2.30GHz
500Gb RAM, 64 CPUs
```
 

### blacklist chimeras

```bash
awk -v "fname1"=${splitdir}/${ext}_${baseoutname}_norm.txt -v "fname2"=${splitdir}/${ext}_${baseoutname}_abnorm.sam -v "fname3"=${splitdir}/${ext}_${baseoutname}_unmapped.sam -f ${juiceDir}/scripts/chimeric_blacklist.awk ${splitdir}/${ext}_${baseoutname}.sam
```

### fragment assignment

```bash
perl ${juiceDir}/scripts/fragment.pl ${splitdir}/${ext}_${baseoutname}_norm.txt ${splitdir}/${ext}_${baseoutname}.frag.txt ${site_file}
```

### sort fragments

sort by chr, frag, strand, and position

```bash
sort --parallel=${threads} \
     --buffer-size=${memory} \
     --temporary-directory=${tmpdir} \
     --key=2,2d --key=6,6d --key=4,4n --key=8,8n --key=1,1n --key=5,5n --key=3,3n ${splitdir}/${ext}_${baseoutname}.frag.txt > ${splitdir}/${ext}_${baseoutname}.sort.txt
```

### merge fragments

```bash
cp ${splitdir}/${ext}_${baseoutname}.sort.txt ${outputdir}/merged_sort.txt
```

### remove duplicates

```bash
awk -v dir=${projectDir} -v groupname=${groupname} -v debugdir=${debugdir} -v juicedir=${juiceDir} -v site=${site} -v genomeID=${genomeID} -v genomePath=${genomePath} -v justexact=0 -f ${juiceDir}/scripts/split_rmdups.awk ${outputdir}/merged_sort.txt > cmds.sh
chmod +x cmds.sh
./cmds.sh
```
If your goal is using the HiC data for geneome assembly, you can stop here and use 
## Step 5. Post analyses stats

```bash
cat ${headfile} > ${outputdir}/inter.txt
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
```

## Step 6. Generate HiC files

```bash
${juiceDir}/scripts/juicer_tools pre \
    -f ${site_file} \
    -s ${outputdir}/inter.txt \
    -g ${outputdir}/inter_hists.m \
    -q 1 ${outputdir}/merged_nodups.txt ${outputdir}/inter.hic ${genomePath}
${juiceDir}/scripts/juicer_tools pre \
    -f ${site_file} \
    -s ${outputdir}/inter_30.txt \
    -g ${outputdir}/inter_30_hists.m \
    -q 30 ${outputdir}/merged_nodups.txt ${outputdir}/inter_30.hic ${genomePath}
```

## Step 7. Other (optional)

Some steps that are in Juicer scripts, if you need them to be run. For genome assemblies, we just need the `merged_nodups.txt` file.

```bash
#collect collisions and dedup
cat $splitdir/*_abnorm.sam > $outputdir/abnormal.sam
cat $splitdir/*_unmapped.sam > $outputdir/unmapped.sam
awk -f ${juiceDir}/scripts/collisions.awk ${outputdir}/abnormal.sam > ${outputdir}/collisions.txt
gawk -v fname=${outputdir}/collisions.txt -f ${juiceDir}/scripts/collisions_dedup_rearrange_cols.awk ${outputdir}/collisions.txt |\
      sort -k3,3n -k4,4n -k10,10n -k11,11n -k17,17n -k18,18n -k24,24n -k25,25n -k31,31n -k32,32n |\
      awk -v name=${outputdir} -f ${juiceDir}/scripts/collisions_dups.awk
# annotating loops
${juiceDir}/scripts/juicer_hiccups.sh \
    -j ${juiceDir}/scripts/juicer_tools \
    -i $outputdir/inter_30.hic \
    -m ${juiceDir}/references/motif \
    -g $genomeID
# annotate contact domains
${juiceDir}/scripts/juicer_arrowhead.sh \
    -j ${juiceDir}/scripts/juicer_tools \
    -i ${outputdir}/inter_30.hic
```