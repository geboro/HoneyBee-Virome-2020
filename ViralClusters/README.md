This folder contains: 
  * The nucleotide sequences of all viral contigs in folder [nucleotideSeqs](./nucleotideSeqs)
  * All protein sequences derived from these viral contigs in file [ViralClusterContigs.faa](./ViralClusterContigs.faa)
  
The header of each protein sequence in the [ViralClusterContigs.faa] file has the following syntax:
  * 1. The Viral Cluster Identifier (e.g. `VC_122_0_`). Unclustered contigs or singletons are marked as `VC_Unclustered`
  * 2. The contig Identifier (e.g. `LD7_45_`), which contains the source of that contig: 
      * `GR7` are contigs assembled from the viral metagenome from Grammont
      * `LD7` are contigs assembled from the viral metagenome from LesDroites
      * `PRO` are contigs identified in bacterial genomes
      * `ISO` are contigs from pure phage isolates 
  * 3. a unique contig identifier from Prokka (e.g. `FMGFAOAF_`)
  * 4. The protein identifier or gene number within that contig (e.g. `00040`)
  * 5. The annotation for that protein sequence assigned by Prokka (mostly `hypothetical protein`)
  
