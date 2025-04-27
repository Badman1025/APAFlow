# APAFlow
The APAFlow pipeline is a comprehensive bioinformatics tool designed to integrate and automate the analysis of 3' UTR alternative polyadenylation (APA) events. It processes DAPARS2 results, performs motif discovery using DREME, and aligns data with TarBase references. This pipeline offers a streamlined and automated approach to studying APA events and their regulatory networks, enhancing efficiency and reproducibility in bioinformatics research.

The APAFlow pipeline is structured into five modules:
* **Process DAPARS2 Results**:Converts DAPARS2 output into a standardized BED file for downstream analysis.
* **Run DREME for Motif Discovery**: Identifies de novo motifs using DREME and aligns them to known motifs via Tomtom for functional annotation.
* **Filter TarBase Interactions**:Extracts experimentally validated 3’UTR miRNA-mRNA interactions from the [TarBase](https://dianalab.e-ce.uth.gr/tarbasev9) database.
* **Integrate Motif and Interaction Data**: Aligns Tomtom motif results with TarBase interactions to identify regulatory overlaps.
* **Generate Interactive Visualizations**: Executes an R script to produce an interactive subnetwork highlighting top-ranked APA regulatory interactions.

## install


## test

