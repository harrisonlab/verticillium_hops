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
for RawData in $(ls qc_dna/paired/*/*/*/*fq.gz | grep '55'); do
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
11055	68.38
11100	61.62
```

# Assembly
Assembly was performed using: Velvet / Abyss / Spades

## Spades Assembly

(Run from the new cluster)

```bash
  for StrainPath in $(ls -d qc_dna/paired/*/*); do
    ProgDir=/projects/oldhome/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades
    Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
    Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
    F_Read=$(ls $StrainPath/F/*.fq.gz)
    R_Read=$(ls $StrainPath/R/*.fq.gz)
    OutDir=assembly/spades/$Organism/${Strain}
    echo $F_Read
    echo $R_Read
    sbatch $ProgDir/slurm_spades.sh $F_Read $R_Read $OutDir correct 10
    # OutDir=assembly/spades/$Organism/${Strain}_2
    # sbatch $ProgDir/slurm_spades.sh $F_Read $R_Read $OutDir correct
  done
```

Quast

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/spades/*/*/filtered_contigs/contigs_min_500bp.fasta); do
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
  # OutDir=assembly/spades/$Organism/$Strain
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

The results of quast were shown using the following commands:

```bash
  for Assembly in $(ls assembly/spades/*/*/filtered_contigs/transposed_report.tsv); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev);
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev);
    # echo;
    # echo $Organism;
    # echo $Strain;
    cat $Assembly | tail -n +2 | sed "s/contigs_min_500bp/${Organism}_${Strain}/g"
  done > assembly/quast_results.txt
```


Contigs were renamed in accordance with ncbi recomendations.

```bash
  ProgDir=$(ls -d /projects/oldhome/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants)
  touch tmp.csv
  for Assembly in $(ls assembly/spades/*/*/filtered_contigs/contigs_min_500bp.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    OutDir=assembly/spades/$Organism/$Strain/ncbi_edits
    mkdir $OutDir
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/contigs_min_500bp_renamed.fasta --coord_file tmp.csv
  done
  rm tmp.csv
```

<!--
A Bioproject and Biosample was made with NCBI genbank for submission of genomes.
Following the creation of these submissions, the .fasta assembly was uploaded
through the submission portal. A note was provided requesting that the assembly
be run through the contamination screen to aid a more detailed resubmission in
future. The returned FCSreport.txt was downloaded from the NCBI webportal and
used to correct the assembly to NCBI standards.

NCBI reports (FCSreport.txt) were manually downloaded to the following loactions:

```bash
  for Assembly in $(ls assembly/spades/*/*/*/contigs_min_500bp_renamed.fasta | grep -v 'ncbi_edits'); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    NCBI_report_dir=genome_submission/$Organism/$Strain/initial_submission
    mkdir -p $NCBI_report_dir
  done
```


These downloaded files were used to correct assemblies:

```bash
for Assembly in $(ls assembly/spades/*/*/ncbi_edits/contigs_min_500bp_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
NCBI_report=$(ls genome_submission/$Organism/$Strain/initial_submission/Contamination*.txt)
OutDir=assembly/spades/$Organism/$Strain/ncbi_edits
mkdir -p $OutDir
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --keep_mitochondria --inp $Assembly --out $OutDir/contigs_min_500bp_renamed.fasta --coord_file $NCBI_report > $OutDir/log.txt
done
```

Quast

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/spades/*/*/ncbi_edits/contigs_min_500bp_renamed.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    OutDir=assembly/spades/$Organism/$Strain/ncbi_edits
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```
-->




# Repeat masking
Repeat masking was performed and used the following programs: Repeatmasker Repeatmodeler

The renamed assembly was used to perform repeatmasking.

```bash
for Assembly in $(ls assembly/spades/*/*/ncbi_edits/contigs_min_500bp_renamed.fasta | grep -v '_2' | grep -v '11055'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism"
echo "$Strain"
OutDir=repeat_masked/$Organism/$Strain/ncbi_edits_repmask
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
qsub $ProgDir/rep_modeling.sh $Assembly $OutDir
qsub $ProgDir/transposonPSI.sh $Assembly $OutDir
done
```

Busco was run to check gene space in assemblies

```bash
  for Assembly in $(ls repeat_masked/*/*/ncbi_edits_repmask/*_contigs_unmasked.fa); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=gene_pred/busco/$Organism/$Strain/assembly
    qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
  done
```

```bash
  for File in $(ls gene_pred/busco/*/*/assembly/*/short_summary_*.txt); do
  Strain=$(echo $File| rev | cut -d '/' -f4 | rev)
  Organism=$(echo $File | rev | cut -d '/' -f5 | rev)
  Complete=$(cat $File | grep "(C)" | cut -f2)
  Single=$(cat $File | grep "(S)" | cut -f2)
  Fragmented=$(cat $File | grep "(F)" | cut -f2)
  Missing=$(cat $File | grep "(M)" | cut -f2)
  Total=$(cat $File | grep "Total" | cut -f2)
  echo -e "$Organism\t$Strain\t$Complete\t$Single\t$Fragmented\t$Missing\t$Total"
  done
```

# Gene prediction


RNAseq data was downloaded for gene annotation.
https://www.ncbi.nlm.nih.gov/bioproject/PRJNA283258
From biosamples:
https://www.ncbi.nlm.nih.gov/sra/SRX1020679[accn]
https://www.ncbi.nlm.nih.gov/sra/SRX1020629[accn]

There is additional RNAseq data from an infection timecourse on hop (https://www.ncbi.nlm.nih.gov/bioproject/PRJEB14243), but this will require more time to map and assemble than we have atm.

```bash
Study=PRJNA283258
Treatment=RECICA91
Rep=1
Accession=SRR2012821
OutDir=raw_rna/$Study/$Treatment/$Rep
mkdir -p $OutDir
fastq-dump --outdir $OutDir --gzip --split-3 $Accession
Rep=2
Accession=SRR2012820
OutDir=raw_rna/$Study/$Treatment/$Rep
mkdir -p $OutDir
fastq-dump --outdir $OutDir --gzip --split-3 $Accession
Rep=3
Accession=SRR2012822
OutDir=raw_rna/$Study/$Treatment/$Rep
mkdir -p $OutDir
fastq-dump --outdir $OutDir --gzip --split-3 $Accession

Study=PRJNA283258
Treatment=T2
Rep=1
Accession=SRR2012772
OutDir=raw_rna/$Study/$Treatment/$Rep
mkdir -p $OutDir
fastq-dump --outdir $OutDir --gzip --split-3 $Accession
Rep=2
Accession=SRR2012773
OutDir=raw_rna/$Study/$Treatment/$Rep
mkdir -p $OutDir
fastq-dump --outdir $OutDir --gzip --split-3 $Accession
Rep=3
Accession=SRR2012794
OutDir=raw_rna/$Study/$Treatment/$Rep
mkdir -p $OutDir
fastq-dump --outdir $OutDir --gzip --split-3 $Accession
```
