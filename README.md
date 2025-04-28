# APAFlow
The APAFlow pipeline is a comprehensive bioinformatics tool designed to integrate and automate the analysis of 3' UTR alternative polyadenylation (APA) events. It processes DAPARS2 results, performs motif discovery using DREME, and aligns data with TarBase references. This pipeline offers a streamlined and automated approach to studying APA events and their regulatory networks, enhancing efficiency and reproducibility in bioinformatics research.

The APAFlow pipeline is structured into five modules:
* **Process DAPARS2 Results**: Converts DAPARS2 output into a standardized BED file for downstream analysis.
* **Run DREME for Motif Discovery**: Identifies de novo motifs using DREME and aligns them to known motifs via Tomtom for functional annotation.
* **Filter TarBase Interactions**: Extracts experimentally validated 3’UTR miRNA-mRNA interactions from the [TarBase](https://dianalab.e-ce.uth.gr/tarbasev9) database.
* **Integrate Motif and Interaction Data**: Aligns Tomtom motif results with TarBase interactions to identify regulatory overlaps.
* **Generate Interactive Visualizations**: Executes an R script to produce an interactive subnetwork highlighting top-ranked APA regulatory interactions.

## install
```
conda env create -f environment.yml
conda activate APAFlow

python APAFlow.py -h
usage: APAFlow.py [-h] -i INPUT_FILE -g GENOME_FASTA -t TARBASE_FILE -m MEME_FILE [-d DISTANCE] [-n TOP_N]

Integrate DAPARS2 result processing, DREME analysis, TarBase filtering, alignment, and subnetwork visualization

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input_file INPUT_FILE
                        Path to the DAPARS/DAPARS2 result file
  -g GENOME_FASTA, --genome_fasta GENOME_FASTA
                        Path to the genome FASTA file
  -t TARBASE_FILE, --tarbase_file TARBASE_FILE
                        Path to the TarBase file
  -m MEME_FILE, --meme_file MEME_FILE
                        Path to the MEME file (e.g., Homo_sapiens_hsa.meme)
  -d DISTANCE, --distance DISTANCE
                        Distance flanking PAS positions (default: 50)
  -n TOP_N, --top_n TOP_N
                        Number of top connections to display (default: 10)
```
## test
- **-i** Testdata(-----------): 3'UTR APA events and distal/proximal PAS sites.  
- **-g** Reference genome: [Ensembl](https://grch37.ensembl.org/index.html/).  
- **-t** Dataset: [TarBasev9](https://dianalab.e-ce.uth.gr/tarbasev9/downloads).  
- **-m** MEME file: miRNA_MEME(----------------).
```
wget https://dianalab.e-ce.uth.gr/tarbasev9/data/Mus_musculus_TarBase-v9.tsv.gz
gunzip Mus_musculus_TarBase-v9.tsv.gz
python APAFlow.py -i all_result_temp.chr2.txt -g Mus_musculus.GRCm39.dna.primary_assembly.fa -t Mus_musculus_TarBase-v9.tsv -m ./miRNA_meme/Mus_musculus_mmu.meme
```














