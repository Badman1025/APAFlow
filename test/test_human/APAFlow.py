import pandas as pd
import argparse
import os
import subprocess

def process_dapars2_results(input_file, distance=50):
    """
    Process DAPARS2 results to generate a BED file.
    """
    df = pd.read_csv(input_file, sep='\t')
    df = df[['Gene', 'Predicted_Proximal_APA']]
    df['Start'] = df['Predicted_Proximal_APA'] - distance
    df['End'] = df['Predicted_Proximal_APA'] + distance
    df = df.drop(columns=['Predicted_Proximal_APA'])
    df[['NM', 'Gene', 'Chr', 'Strand']] = df['Gene'].str.split('|', expand=True)
    df = df.drop(columns=['NM'])
    df['Chr'] = df['Chr'].str.replace('chr', '')
    df = df[['Chr', 'Start', 'End', 'Gene', 'Strand']]
    output_bed_file = "3UTRAPA.bed"
    df.to_csv(output_bed_file, sep='\t', index=False, header=False)
    print(f"Processing completed. Results saved to {output_bed_file}")
    return output_bed_file

def run_dreme(genome_fasta, meme_file):
    """
    Run DREME tool.
    """
    print("Running DREME tool...")
    bed_file = "3UTRAPA.bed"
    fasta_file = "3UTRAPA.fa"
    subprocess.run(f"bedtools getfasta -fi {genome_fasta} -bed {bed_file} -fo {fasta_file}", shell=True)
    dreme_output_dir = "dreme_3UTRrna"
    subprocess.run(f"rm -rf {dreme_output_dir}", shell=True)
    subprocess.run(f"dreme -o {dreme_output_dir} -p {fasta_file} -rna", shell=True)
    
    # Generate target.meme in the current directory
    subprocess.run(f"meme2meme {dreme_output_dir}/dreme.txt > target.meme", shell=True)
    
    # Remove existing tomtom_out directory if it exists
    tomtom_output_dir = "tomtom_out"
    subprocess.run(f"rm -rf {tomtom_output_dir}", shell=True)
    
    # Run tomtom and let it create the tomtom_out directory automatically
    subprocess.run(f"tomtom target.meme {meme_file} -o {tomtom_output_dir}", shell=True)
    print("DREME tool execution completed")
    
    tomtom_file = f"{tomtom_output_dir}/tomtom.tsv"
    if not os.path.exists(tomtom_file):
        print("The motif was not found in this region.")
        return None
    return tomtom_file

def filter_tarbase(tarbase_file):
    """
    Filter TarBase file to extract 3UTR data.
    """
    df = pd.read_csv(tarbase_file, sep='\t')
    first_row = df.iloc[0:1]
    filtered_df = df[df.iloc[:, 5] == '3UTR']
    result_df = pd.concat([first_row, filtered_df], ignore_index=True)
    tarbase_3utr_file = "TarBase_3UTR.tsv"
    result_df.to_csv(tarbase_3utr_file, sep='\t', index=False)
    print(f"Filtering completed. Results saved to {tarbase_3utr_file}")
    return tarbase_3utr_file

def align_tomtom_tarbase(tomtom_file, tarbase_3utr_file):
    """
    Align Tomtom and TarBase files.
    """
    if tomtom_file is None:
        print("Skipping alignment due to missing Tomtom results.")
        return

    tomtom_df = pd.read_csv(tomtom_file, sep='\t')
    filtered_tomtom_df = tomtom_df[tomtom_df['p-value'] < 0.05]
    tarbase_df = pd.read_csv(tarbase_3utr_file, sep='\t')
    merged_df = pd.merge(filtered_tomtom_df[['Target_ID']], tarbase_df, left_on='Target_ID', right_on='mirna_name')
    bed_df = pd.read_csv("3UTRAPA.bed", sep='\t', header=None, names=['chromosome', 'start', 'end', 'gene_name', 'strand'])
    final_merged_df = pd.merge(merged_df, bed_df, on='gene_name')
    final_output_file = "TarBase_3UTRAPA.tsv"
    final_merged_df.to_csv(final_output_file, sep='\t', index=False)
    print(f"Alignment completed. Results saved to {final_output_file}")
    return final_output_file

def run_r_script(top_n):
    """
    Run the R script for subnetwork visualization.
    """
    print("Running R script for subnetwork visualization...")
    r_script = f"""
    library(igraph)
    library(edgebundleR)
    library(dplyr)
    library(purrr)
    suppressPackageStartupMessages(library(widgetframe))
    library(htmlwidgets)

    parameters <- read.delim("./TarBase_3UTRAPA.tsv", sep = "\\t", header = TRUE)
    parameters$microt_score[is.na(parameters$microt_score)] <- 0
    parameters <- parameters[parameters$microt_score >= 1, ]
    relationship <- data.frame(
      from = parameters$mirna_name,
      to = parameters$gene_name
    )

    # Display top connections
    top_to <- head(sort(table(relationship$to), decreasing = TRUE), n = {top_n})
    top_from <- head(sort(table(relationship$from), decreasing = TRUE), n = {top_n})
    print(head(sort(table(relationship$to), decreasing = TRUE), n = {top_n}))
    print(head(sort(table(relationship$from), decreasing = TRUE), n = {top_n}))

    # Extract top nodes
    top_nodes <- c(names(top_to), names(top_from))

    # Filter relationship to include only top nodes
    filtered_relationship <- relationship %>%
      filter(from %in% top_nodes | to %in% top_nodes)

    # Create graph from filtered data
    graph <- graph_from_data_frame(filtered_relationship, directed = FALSE)

    # Visualize the graph
    widget <- edgebundle(graph, tension = 0.45, cutoff = 0.1, width = NULL, fontsize = 10,
                         padding = 150, nodesize = c(3, 10), directed = FALSE)
    htmlwidgets::saveWidget(widget, "network_visualization.html")

    # Save filtered parameters
    write.csv(parameters, "parameters_microt_score=1.csv", row.names = FALSE)
    """
    with open("subnetwork_visualization.R", "w") as f:
        f.write(r_script)
    subprocess.run("Rscript subnetwork_visualization.R", shell=True)
    print("R script execution completed. Visualization saved to network_visualization.html")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Integrate DAPARS2 result processing, DREME analysis, TarBase filtering, alignment, and subnetwork visualization")
    parser.add_argument("-i", "--input_file", required=True, help="Path to the DAPARS/DAPARS2 result file")
    parser.add_argument("-g", "--genome_fasta", required=True, help="Path to the genome FASTA file")
    parser.add_argument("-t", "--tarbase_file", required=True, help="Path to the TarBase file")
    parser.add_argument("-m", "--meme_file", required=True, help="Path to the MEME file (e.g., Homo_sapiens_hsa.meme)")
    parser.add_argument("-d", "--distance", type=int, default=50, help="Distance flanking PAS positions (default: 50)")
    parser.add_argument("-n", "--top_n", type=int, default=10, help="Number of top connections to display (default: 10)")
    args = parser.parse_args()

    # Block 1: Process DAPARS2 results
    bed_file = process_dapars2_results(args.input_file, args.distance)

    # Block 2: Run DREME tool
    tomtom_file = run_dreme(args.genome_fasta, args.meme_file)

    # Block 3: Filter TarBase file
    tarbase_3utr_file = filter_tarbase(args.tarbase_file)

    # Block 4: Align Tomtom and TarBase files
    tarbase_3utrapa_file = align_tomtom_tarbase(tomtom_file, tarbase_3utr_file)

    # Block 5: Run R script for subnetwork visualization
    if tarbase_3utrapa_file:
        run_r_script(args.top_n)