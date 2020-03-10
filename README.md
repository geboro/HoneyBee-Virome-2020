# HoneyBee-Virome-2020
Pipeline for the analysis of Bonilla-Rosso et al. 2020. "Honey bees harbor a diverse gut virome engagin in nested strain-level interactions with the microbiota" PNAS. [doi:xxxx]

### Assembly
The metaviromes were assembled with **SPAdes** (v3.13.0 downloaded October 2018) [http://cab.spbu.ru/software/spades/] as follows:
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
And evaluated with **QUAST** (v5.0.1 downloaded October 2018) [https://github.com/ablab/quast] as follows:
```bash
quast.py 
  -o ./metaquast01 
  -R <fasta reference genomes> 
  --contig-thresholds 1000,2000,2500,3000,5000,7500,10000,15000,20000,25000,30000,35000,40000,45000,50000,55000,60000,70000,80000,90000,100000 
  -t 25 
  --unique-mapping 
  Grammont01/contigs.fasta  
  LeDroites01/contigs.fasta 
```
### Detection of Viral Contigs
We used **VirSorter** (v1.0.5 downloaded in March 2019) [https://github.com/simroux/VirSorter] as follows:
```bash
wrapper_phage_contigs_sorter_iPlant.pl 
  -f Grammont_assembly.fasta 
  -db 1 
  --wdir /home/gbonilla/Documents/phages/spades/NEW/04_virsorter/Grammont_virsorter
  --virome 
  --ncpu 30 
  --data-dir ./virsorter-data
```

### Detection of Viral Contigs
We used **VirSorter** (v1.0.5 downloaded in March 2019) [https://github.com/simroux/VirSorter] as follows:
```bash
wrapper_phage_contigs_sorter_iPlant.pl 
  -f Grammont_assembly.fasta 
  -db 1 
  --wdir /home/gbonilla/Documents/phages/spades/NEW/04_virsorter/Grammont_virsorter
  --virome 
  --ncpu 30 
  --data-dir ./virsorter-data
```
