Verticillium nonalfalfae ex. hop
====================

Commands used during analysis of the Verticillium nonalfalfae ex. hop genomes. Note - all this work was performed in the directory: /projects/verticillium_hops

The following is a summary of the work presented in this Readme:
Data organisation:
  * Preparing data  
Draft Genome assembly
  * Data qc
  * Genome assembly
  * Repeatmasking
  * Gene prediction
  * Functional annotation
Genome analysis
  *


#Data organisation

Data was copied from the raw_data repository to a local directory for assembly
and annotation.

```bash
    mkdir -p /projects/verticillium_hops
  	cd /projects/verticillium_hops
    RawDat=$(ls -d /archives/2019_niabemr_miseq/RAW/191010_M04465_0101_000000000-C6N3P/Data/Intensities/BaseCalls)
    # Isolate 11043
  	Species=V.nonalfalfae
  	Strain=11043
    OutDir=raw_dna/paired/$Species/$Strain
    mkdir -p $OutDir/F
    mkdir -p $OutDir/R
    cp -s $RawDat/${Strain}_*_R1_001.fastq.gz $PWD/$OutDir/F/.
    cp -s $RawDat/${Strain}_*_R2_001.fastq.gz $PWD/$OutDir/R/.
    # Isolate 11055
    Species=V.nonalfalfae
    Strain=11055
    OutDir=raw_dna/paired/$Species/$Strain
    mkdir -p $OutDir/F
    mkdir -p $OutDir/R
    cp -s $RawDat/${Strain}_*_R1_001.fastq.gz $PWD/$OutDir/F/.
    cp -s $RawDat/${Strain}_*_R2_001.fastq.gz $PWD/$OutDir/R/.
    # Isolate 11100
    Species=V.nonalfalfae
    Strain=11100
    OutDir=raw_dna/paired/$Species/$Strain
    mkdir -p $OutDir/F
    mkdir -p $OutDir/R
    cp -s $RawDat/${Strain}_*_R1_001.fastq.gz $PWD/$OutDir/F/.
    cp -s $RawDat/${Strain}_*_R2_001.fastq.gz $PWD/$OutDir/R/.
```


#Data qc

programs: fastqc fastq-mcf kmc

Data quality was visualised using fastqc:

```bash
  for RawData in $(ls raw_dna/paired/*/*/*/*.fastq.gz); do
    echo $RawData;
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
    qsub $ProgDir/run_fastqc.sh $RawData;
  done
```

Trimming was performed on data to trim adapters from sequences and remove poor quality data.
This was done with fastq-mcf

```bash
for StrainPath in $(ls -d raw_dna/paired/*/*); do
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
	IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
	ReadsF=$(ls $StrainPath/F/*.fastq*)
	ReadsR=$(ls $StrainPath/R/*.fastq*)
	echo $ReadsF
	echo $ReadsR
	qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters DNA
done
```
Data quality was visualised once again following trimming:

```bash
  for RawData in $(ls qc_dna/paired/*/*/*/*.fq.gz); do
    echo $RawData;
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
    qsub $ProgDir/run_fastqc.sh $RawData;
  done
```


Sequencing coverage was estimated:

```bash
for RawData in $(ls qc_dna/paired/*/*/*/*fq.gz); do
echo $RawData;
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
qsub $ProgDir/run_fastqc.sh $RawData
GenomeSz=35
OutDir=$(dirname $RawData)
qsub $ProgDir/sub_count_nuc.sh $GenomeSz $RawData $OutDir
done
```


Find predicted coverage for these isolates:

```bash
for StrainDir in $(ls -d qc_dna/paired/*/*); do
Strain=$(basename $StrainDir)
printf "$Strain\t"
for File in $(ls qc_dna/paired/*/"$Strain"/*/*.txt); do
echo $(basename $File);
cat $File | tail -n1 | rev | cut -f2 -d ' ' | rev;
done | grep -v '.txt' | awk '{ SUM += $1} END { print SUM }'
done
```

```
11043	66.48
11100	61.62
```

# Assembly
Assembly was performed using: Velvet / Abyss / Spades

## Spades Assembly

```bash
  for StrainPath in $(ls -d qc_dna/paired/*/* | tail -n +2); do
    ProgDir=/projects/oldhome/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades
    Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
    Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
    F_Read=$(ls $StrainPath/F/*.fq.gz)
    R_Read=$(ls $StrainPath/R/*.fq.gz)
    OutDir=assembly/spades/$Organism/$Strain
    echo $F_Read
    echo $R_Read
    sbatch $ProgDir/slurm_spades.sh $F_Read $R_Read $OutDir correct 10
  done
```

Quast

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
    for Assembly in $(ls assembly/spades/*/*/ncbi_edits/contigs_min_500bp.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    OutDir=assembly/spades/$Organism/$Strain/ncbi_edits
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

The results of quast were shown using the following commands:

```bash
  for Assembly in $(ls assembly/spades/*/*/ncbi_edits/report.txt); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev);
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev);
    echo;
    echo $Organism;
    echo $Strain;
    cat $Assembly;
  done > assembly/quast_results.txt
```
