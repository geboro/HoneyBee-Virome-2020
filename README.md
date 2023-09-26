# HoneyBee-Virome-2020
Pipeline for the analysis of Bonilla-Rosso et al. 2020. "Honey bees harbor a diverse gut virome engaging in nested strain-level interactions with the microbiota" PNAS [doi:10.1073/pnas.2000228117](https://www.pnas.org/doi/10.1073/pnas.2000228117).

## Data Deposition
Sequencing data produced during this project is deposited at DDBJ/ENA/GenBank under BioProject [PRJNA599270](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA599270). Raw reads are deposited in the SRA under accesions [SRX7509864](https://www.ncbi.nlm.nih.gov/sra/SRX7509864[accn]) for Grammont, and [SRX7509865](https://www.ncbi.nlm.nih.gov/sra/SRX7509865[accn]) for LesDroites. The assembled viral contigs are deposited as a Whole Genome shotgun project under accession JAAOBB000000000. The genomes of the pure phage isolates are deposited 
in GenBank under accesion numbers [MT006233-MT006240](https://www.ncbi.nlm.nih.gov/nuccore?term=MT006233%5Baccn%5D%3AMT006240%5Baccn%5D&cmd=DetailsSearch). 


### Assembly
The metaviromes were assembled with **SPAdes** (v3.13.0) [http://cab.spbu.ru/software/spades/] as follows:
```bash
spades.py 
  -o ./virome 
  --pe1-1 Virome_R1_trim_pair.fastq 
  --pe1-2 Virome_R2_trim_pair.fastq 
  -k 21,33,55,77,99,127 
  --meta 
  -m 120 
  -t 40
  --phred-offset 33
```
And evaluated with **QUAST** (v5.0.1) [https://github.com/ablab/quast] as follows:
```bash
quast.py 
  -o ./metaquast01 
  -R <fasta reference genomes> 
  --contig-thresholds 1000,2000,2500,3000,5000,7500,10000,15000,20000,25000,30000,35000,40000,45000,50000,55000,60000,70000,80000,90000,100000 
  -t 25 
  --unique-mapping 
  Virome/contigs.fasta  
  
```
### Detection of Viral Contigs
We used **VirSorter** (v1.0.5) [https://github.com/simroux/VirSorter] as follows:
```bash
wrapper_phage_contigs_sorter_iPlant.pl 
  -f Virome_assembly.fasta 
  -db 1 
  --wdir /home/gbonilla/Documents/phages/spades/NEW/04_virsorter/Virome_virsorter
  --virome 
  --ncpu 30 
  --data-dir ./virsorter-data
```

### Read Mapping
We used **BBMap** (v37.32) [https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/] to map virome reads to assembled contigs as follows:
```bash
~/bbmap/bbsplit.sh 
  build=1 ref=./virome/infile
~/bbmap/bbsplit.sh 
  build=1
  in=Virome_R1.fastq
  in2=Virome_R2.fastq
  ambiguous2=best 
  basename=%_BBsplit.sam 
  refstats=VR_refstats.OUT 
  nzo=f
```
And **samtools** (v1.9 using htslib 1.9) [https://github.com/samtools/samtools] to generate coverage depth profiles
```bash
for i in *sam; do
  samtools view -bh $i | samtools sort > ${i%.sam}_sorted.bam
  samtools index ${i%.sam}_sorted.bam
  samtools depth -a ${i%.sam}_sorted.bam > ${i%_bbsplitted.sam}_DEPTH.txt
  samtools idxstats ${i%.sam}_sorted.bam > ${i%_bbsplitted.sam}_IDXSTATS.txt
done
```
### Clustering 
We used **vConTACT** (v2.9.10)[https://bitbucket.org/MAVERICLab/vcontact2/src/master/] to cluster viral contigs as follows:
```bash
vcontact 
  -f 
  --raw-proteins virome_FINAL.prt 
  --rel-mode 'Diamond' 
  --proteins-fp prot_2_ctg_FINAL.csv 
  --db 'ProkaryoticViralRefSeq85-Merged' 
  --pcs-mode MCL 
  --vcs-mode ClusterONE 
  --c1-bin /usr/local/bin/cluster_one-1.0.jar 
  --output-dir vcontact_FINAL 
  --vc-penalty 1
```
The similarity network, analysed with **Cytoscape** (v3.7.2), is available as a cytoscape object in [vConTACT_similarity_network.cys](./vConTACT_similarity_network.cys)

### ViralCluster Analysis
We calculated the Average Nucleotide and Aminoacid Identity for each ViralCluster with **enveomics** (v.1.1.4) [https://github.com/lmrodriguezr/enveomics], and analysed as described in [ANIanalyses](ANIanalyses).

### Read Recruitment to ViralClusters
We used **BBMap** (v37.32) [https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/] to map virome reads against contigs grouped in Viral Clusters...
```bash
bbsplit.sh 
  build=175
  in=~/Documents/phages/spades/NEW/02_decontaminate/Grammon2017_R1_decont.fastq
  in2=~/Documents/phages/spades/NEW/02_decontaminate/Grammon2017_R2_decont.fastq
  ambiguous2=all
  minid=0.9
  basename=%_BBsplit.sam 
  bs=GR_shellscriptfile.OUT 
  refstats=GR_refstats.OUT 
  nzo=f
  threads=40

bbsplit.sh 
  build=175
  in=~/Documents/phages/spades/NEW/02_decontaminate/LeDroite2018_R1_decont.fastq
  in2=~/Documents/phages/spades/NEW/02_decontaminate/LeDroite2018_R2_decont.fastq
  ambiguous2=all
  minid=0.9
  basename=%_BBsplit.sam 
  bs=LD_shellscriptfile.OUT 
  refstats=LD_refstats.OUT 
  nzo=f
  threads=40
```
And then caclulated skewness with [ctgEvenness.py](https://github.com/geboro/HoneyBee-Virome-2020/blob/master/ctgEvenness.py):
```bash
for i in *sam; do
  samtools view -bh $i | samtools sort > ${i%.sam}_sorted.bam
  samtools index ${i%.sam}_sorted.bam
  samtools depth -a ${i%.sam}_sorted.bam > ${i%_bbsplitted.sam}_DEPTH.txt
  samtools idxstats ${i%.sam}_sorted.bam > ${i%_bbsplitted.sam}_IDXSTATS.txt
  echo "contig  ctgLen  totReads(CoverageDensity)   coverdBases percGenomeCov   covMedian   obsEve  equitability    covCVar covKurto    covSkew"
  python ctgEvenness.py $i >> ${i}.ctgStats
done
cat *ctgStats > 0Mapping_Summary_GR.tsv

#and then calculate the number of bases within each vcluster:
cd ~/Documents/phages/spades/NEW/09_vcmaps/infiles
for i in *fasta; do 
  echo -n "${i} "
  grep -v ">" ${i} | wc | awk '{print $3-$1}'
done
```

### Mapping reads from Bacterial Metagenomes
We used **BBMap** (v37.32) [https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/] to map reads from 73 bacterial metagenomes against contigs grouped in Viral Clusters...

```bash
for i in ~/KEmetagenomes/*R1.fastq.gz; do
  long=${i%_unmapped_R[12].fastq.gz}; 
  name=${long#*KEmetagenomes/};
  mkdir ./mapping_files/${name};
  bbsplit.sh 
    build=175 
    in=${i} 
    in2=${i%1.fastq.gz}2.fastq.gz 
    ambiguous2=all 
    basename=%_${name}_BBsplit.sam 
    refstats=${name}_refstats.OUT 
    nzo=f 
    threads=40;  
  cat *OUT | awk '{sum+=$6} END { print "AVG_unamb_READS = ",sum/NR}' > ./mapping_files/0_${name}_stats
  cat *OUT | awk '{sum+=$7} END { print "AVG_amb_READS = ",sum/NR}' >> ./mapping_files/0_${name}_stats
  cat *OUT | awk '{sum+=$7} END { print "TOT_amb_READS = ",sum}' >> ./mapping_files/0_${name}_stats
  cat *OUT | awk '{sum+=$6} END { print "TOT_unamb_READS = ",sum}' >> ./mapping_files/0_${name}_stats
  mv *sam ./mapping_files/${name};
  mv *OUT ./mapping_files/${name};
done
```
