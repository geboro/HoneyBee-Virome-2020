# HoneyBee-Virome-2020
Pipeline for the analysis of Bonilla-Rosso et al. 2020. "Honey bees harbor a diverse gut virome engagin in nested strain-level interactions with the microbiota" PNAS. [doi:xxxx]

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
### ViralCluster Analysis
We calculated the Average Nucleotide and Aminoacid Identity for each ViralCluster with **enveomics** (v.1.1.4) [https://github.com/lmrodriguezr/enveomics], and edited
