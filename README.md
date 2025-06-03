# APAFlow
The APAFlow pipeline is a comprehensive bioinformatics tool designed to integrate and automate the analysis of 3' UTR alternative polyadenylation (APA) events. It processes Dapars results, performs motif discovery using DREME, and aligns data with TarBase references. This pipeline offers a streamlined and automated approach to studying APA events and their regulatory networks, enhancing efficiency and reproducibility in bioinformatics research.
![image](https://github.com/user-attachments/assets/746220f7-3ca5-45be-8136-68436f82c028)


The APAFlow pipeline is structured into five modules:
* **Process Dapars Results**: Converts Dapars/Dapars2/Dapars_LR output into a standardized BED file for downstream analysis.
* **Run DREME for Motif Discovery**: Identifies de novo motifs using DREME and aligns them to known motifs via Tomtom for functional annotation.
* **Filter TarBase Interactions**: Extracts experimentally validated 3’UTR miRNA-mRNA interactions from the [TarBase](https://dianalab.e-ce.uth.gr/tarbasev9) database.
* **Integrate Motif and Interaction Data**: Aligns Tomtom motif results with TarBase interactions to identify regulatory overlaps.
* **Generate Interactive Visualizations**: Executes an R script to produce an interactive subnetwork highlighting top-ranked APA regulatory interactions.

## Install
```
conda env create -f environment.yml
conda activate APAFlow
R
install.packages("edgebundleR")

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
## Test
- **`-i`** [Testdata](https://github.com/Badman1025/APAFlow/tree/main/test): 3'UTR APA events and distal/proximal PAS.  
- **`-g`** Reference genome: [Ensembl](https://grch37.ensembl.org/index.html/).  
- **`-t`** Dataset: [TarBasev9](https://dianalab.e-ce.uth.gr/tarbasev9/downloads).  
- **`-m`** MEME file: [miRNA_MEME](https://github.com/Badman1025/APAFlow/tree/main/miRNA_meme).
```
wget https://dianalab.e-ce.uth.gr/tarbasev9/data/Mus_musculus_TarBase-v9.tsv.gz
wget https://dianalab.e-ce.uth.gr/tarbasev9/data/Homo_sapiens_TarBase-v9.tsv.gz

python APAFlow.py -i all_result_temp.chr2.txt -g Mus_musculus.GRCm39.dna.primary_assembly.fa -t Mus_musculus_TarBase-v9.tsv -m ./miRNA_meme/Mus_musculus_mmu.meme
python APAFlow.py -i all_result_temp.chr1.txt -g Homo_sapiens.GRCh38.dna.primary_assembly.fa -t Homo_sapiens_TarBase-v9.tsv -m ./miRNA_meme/Homo_sapiens_hsa.meme -d 50 -n 1
```
## Result
![image](https://github.com/user-attachments/assets/32f4c214-dae1-49a8-a10f-1ed437339994)

## Custom Operation
After running APAFlow.py, the file `TarBase_3UTRAPA.tsv` will be generated. If you need to focus on a specific type of tissue/cell, you can manually filter the table. The filtered file, `TarBase_3UTRAPA_selected.tsv`, can then be processed using [Subnetwork_Visualization_selected.R](https://github.com/Badman1025/APAFlow/blob/main/subnetwork_visualization.R). This script allows filtering based on the `microt_score` parameter and visualizes the entire/top n subnetwork graphs.
```
Rscript Subnetwork_Visualization_selected.R
```
## Citation
1. article
2. Giorgos Skoufos, Panos Kakoulidis, Spyros Tastsoglou, Elissavet Zacharopoulou, Vasiliki Kotsira, Marios Miliotis, Galatea Mavromati, Dimitris Grigoriadis, Maria Zioga, Angeliki Velli, Ioanna Koutou, Dimitra Karagkouni, Steve Stavropoulos, Filippos S Kardaras, Anna Lifousi, Eustathia Vavalou, Armen Ovsepian, Anargyros Skoulakis, Sotiris K Tasoulis, Spiros V Georgakopoulos, Vassilis P Plagianakos, Artemis G Hatzigeorgiou, TarBase-v9.0 extends experimentally supported miRNA–gene interactions to cell-types and virally encoded miRNAs, Nucleic Acids Research, 2023, DOI:  https://doi.org/10.1093/nar/gkad1071
3. Timothy L. Bailey, DREME: motif discovery in transcription factor ChIP-seq data, Bioinformatics, Volume 27, Issue 12, June 2011, Pages 1653–1659, https://doi.org/10.1093/bioinformatics/btr261


