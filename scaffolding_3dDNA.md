# Scaffolding Draft genome

Using the HiC data processed with Juicer, you can scaffold a draft, contig-level genome assembly. You will need:

1. Contig assembly
2. Processed HiC data (`merged_nodups.txt`)


## Enviornmnet setup

```bash
mkdir -p scaffold_hiC
cd scaffold_hiC
git clone git@github.com:HuffordLab/Juicer_PanAnd.git
mv Juicer_PanAnd juicerdir
cd juicerdir
conda env create -f 3ddna.yml
conda activate 3DDNA
cd ..
git clone https://github.com/aidenlab/3d-dna.git
cd 3d-dna
export PATH=$(pwd):$PATH
chmod +x *.sh
```


## Required files  

You can either copy/move/softlink your files to this directory or run them in the same directory where you ran `juicer`, it is up to you. Here, we softlink files required to run 3D-DNA.

```bash
# be sure to use the same genome you used for juicer
ln -s /path/to/your/contigs.fasta ./scaffold_hiC/
ln -s /path/to/your/merged_nodups.txt ./scaffold_hiC/
```


## Run 3D-DNA scaffolder

There are many settings for scaffolding, be sure to check with `--help` and use the ones that are relavant for your case.

```bash
cd scaffold_hiC
run-asm-pipeline.sh \
    --editor-repeat-coverage 20 \
    --rounds 5 \
    --mode haploid \
    contigs.fasta \
    merged_nodups.txt
```
