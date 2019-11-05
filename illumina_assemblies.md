Verticillium nonalfalfae ex. hop
====================

Commands used during analysis of the Verticillium nonalfalfae ex. hop genomes. Note - all this work was performed in the directory: /projects/verticillium_hops

The following is a summary of the work presented in this readme:
Data organisation:
  * Preparing data  
  * Draft Genome assembly
  * Data qc
  * Genome assembly
  * Repeatmasking
  * Gene prediction
  * Functional annotation
  * Genome analysis


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

```
V.nonalfalfae_11043	776	654	33594125	33506343	776	520303	33594125	55.47	104159	61674	93	197	3.65
V.nonalfalfae_11055	68	48	4757506	4743134	68	726754	4757506	65.67	330854	207683	6	10	0.00
V.nonalfalfae_11100	698	590	33233147	33155996	698	392287	33233147	55.85	117696	63755	91	185	3.5
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


The TransposonPSI masked bases were used to mask additional bases from the
repeatmasker / repeatmodeller softmasked and hardmasked files.

```bash

for File in $(ls repeat_masked/*/*/*/*_contigs_softmasked.fa); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
echo "Number of masked bases:"
cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
done
# The number of N's in hardmasked sequence are not counted as some may be present within the assembly and were therefore not repeatmasked.
for File in $(ls repeat_masked/*/*/*/*_contigs_hardmasked.fa); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_hardmasked.fa/_contigs_hardmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -fi $File -bed $TPSI -fo $OutFile
done
```

```
repeat_masked/V.nonalfalfae/11043/ncbi_edits_repmask/11043_contigs_softmasked_repeatmasker_TPSI_appended.fa
Number of masked bases:
1263803
repeat_masked/V.nonalfalfae/11100/ncbi_edits_repmask/11100_contigs_softmasked_repeatmasker_TPSI_appended.fa
Number of masked bases:
772617
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

```
V.nonalfalfae	11043	3663	3656	29	33	3725
V.nonalfalfae	11100	3673	3666	23	29	3725
```

# Gene prediction


RNAseq data was downloaded for gene annotation.
https://www.ncbi.nlm.nih.gov/bioproject/PRJNA283258
From biosamples:
https://www.ncbi.nlm.nih.gov/sra/SRX1020679[accn]
https://www.ncbi.nlm.nih.gov/sra/SRX1020629[accn]

There is additional RNAseq data from an infection timecourse on hop (https://www.ncbi.nlm.nih.gov/bioproject/PRJEB14243), but this will require more time to map and assemble than we have atm.

<!--
```bash
TmpDir=/projects/public_sra
mkdir -p $TmpDir
chmod a+rw $TmpDir
mkdir $HOME/.ncbi
echo "/repository/user/main/public/root = \"$TmpDir\"" > $HOME/.ncbi/user-settings.mkfg
```
-->

```bash
Study=PRJNA283258
Treatment=RECICA91
Rep=1
Accession=SRR2012821
OutDir=raw_rna/$Study/$Treatment/$Rep
mkdir -p $OutDir
# fastq-dump --outdir $OutDir --gzip --split-3 $Accession
fasterq-dump --outdir $OutDir --split-3 -verbose $Accession
Rep=2
Accession=SRR2012820
OutDir=raw_rna/$Study/$Treatment/$Rep
mkdir -p $OutDir
fasterq-dump --outdir $OutDir --split-3 -verbose $Accession
Rep=3
Accession=SRR2012822
OutDir=raw_rna/$Study/$Treatment/$Rep
mkdir -p $OutDir
fasterq-dump --outdir $OutDir --split-3 -verbose $Accession

Study=PRJNA283258
Treatment=T2
Rep=1
Accession=SRR2012772
OutDir=raw_rna/$Study/$Treatment/$Rep
mkdir -p $OutDir
fasterq-dump --outdir $OutDir --split-3 -verbose $Accession
Rep=2
Accession=SRR2012773
OutDir=raw_rna/$Study/$Treatment/$Rep
mkdir -p $OutDir
fasterq-dump --outdir $OutDir --split-3 -verbose $Accession
Rep=3
Accession=SRR2012794
OutDir=raw_rna/$Study/$Treatment/$Rep
mkdir -p $OutDir
fasterq-dump --outdir $OutDir --split-3 -verbose $Accession
```

The reference assembly and gene models were downloaded for V. nonalfalfae genome VnAa140:
https://www.ncbi.nlm.nih.gov/genome/genomes/43074?

These were uploaded from my local computer into the following directory:
```bash
mkdir -p assembly/external_groups/V.nonalfalfae/VnAa140/ncbi
```

```bash
scp /Users/armita/Downloads/vert_hops/* cluster:/projects/verticillium_hops/assembly/external_groups/V.nonalfalfae/VnAa140/ncbi/.
```


### QC of RNAseq data:

```bash
for StrainPath in $(ls -d raw_rna/*/*/*); do
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
	IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
	ReadsF=$(ls $StrainPath/*_1.fastq)
	ReadsR=$(ls $StrainPath/*_2.fastq)
	echo $ReadsF
	echo $ReadsR
	qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters RNA
done
```

## Alignment

Then Rnaseq data was aligned to each genome assembly:

```bash
for Assembly in $(ls repeat_masked/*/*/ncbi_edits_repmask/*_contigs_unmasked.fa); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
for FileF in $(ls qc_rna/raw_rna/*/*/F/*_trim.fq.gz); do
  FileR=$(echo $FileF | sed 's&/F/&/R/&g' | sed 's/_1_trim/_2_trim/g')
  Study=$(echo $FileF | rev | cut -f3 -d '/' | rev)
  Timepoint=$(echo $FileF | rev | cut -f3 -d '/' | rev)
OutDir=alignment/$Organism/$Strain/$Study/$Timepoint
# Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
# while [ $Jobs -gt 1 ]; do
# sleep 1m
# printf "."
# Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
# done
# printf "\n"
# printf "\n"
echo "$Organism - $Strain - $Study - $Timepoint"
# echo $FileF
# echo $FileR
# Prefix=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_star.sh $Assembly $FileF $FileR $OutDir
echo "$Strain\t$Timepoint" >> alignment.log
done
done
```

```bash
mkdir qc_rna/concat
cat qc_rna/raw_rna/*/*/F/*.fq.gz > qc_rna/concat/concat_1.fq.gz
cat qc_rna/raw_rna/*/*/R/*.fq.gz > qc_rna/concat/concat_2.fq.gz

for Assembly in $(ls repeat_masked/*/*/ncbi_edits_repmask/*_contigs_unmasked.fa | head -n1); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
FileF=qc_rna/concat/concat_1.fq.gz
FileR=qc_rna/concat/concat_2.fq.gz
OutDir=alignment/$Organism/$Strain/concat
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_star.sh $Assembly $FileF $FileR $OutDir
done
```

## FunGap

Gene prediction was run on the new cluster from a screen session with a ssh connection to a worker node:

```bash
screen -a

gunzip qc_rna/concat/concat_*.fq.gz
conda activate fungap

for Assembly in $(ls repeat_masked/*/*/ncbi_edits_repmask/*_contigs_unmasked.fa); do
# for Assembly in $(ls repeat_masked/*/*/ncbi_edits_repmask/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
FileF=$(ls qc_rna/concat/concat_1.fq)
FileR=$(ls qc_rna/concat/concat_2.fq)
RefProteome=$(ls assembly/external_groups/V.nonalfalfae/VnAa140/ncbi/GCA_003724135.1_ASM372413v1_protein.faa
)
ReadF=$(ls qc_rna/concat/concat_1.fq)
ReadR=$(ls qc_rna/concat/concat_2.fq)
BamAlignment=$(ls alignment/$Organism/$Strain/concat/star_aligmentAligned.sortedByCoord.out.bam)
Prefix=${Organism}_${Strain}
OutDir=gene_pred/fungap/$Organism/${Strain}
# OutDir=gene_pred/fungap/$Organism/${Strain}_softmasked
ProgDir=/projects/oldhome/armita/git_repos/emr_repos/tools/gene_prediction/fungap
sbatch $ProgDir/slurm_fungap.sh $Assembly $RefProteome $ReadF $ReadR $BamAlignment $Prefix $OutDir
done
```
