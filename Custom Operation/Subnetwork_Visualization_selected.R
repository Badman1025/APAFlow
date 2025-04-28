
    library(igraph)
    library(edgebundleR)
    library(dplyr)
    library(purrr)
    suppressPackageStartupMessages(library(widgetframe))
    library(htmlwidgets)

    parameters <- read.delim("./TarBase_3UTRAPA_selected.tsv", sep = "\t", header = TRUE)
    parameters$microt_score[is.na(parameters$microt_score)] <- 0
    parameters <- parameters[parameters$microt_score >= 1, ]
    relationship <- data.frame(
      from = parameters$mirna_name,
      to = parameters$gene_name
    )

    # Display top connections
    top_to <- head(sort(table(relationship$to), decreasing = TRUE), n = 5)
    top_from <- head(sort(table(relationship$from), decreasing = TRUE), n = 5)
    print(head(sort(table(relationship$to), decreasing = TRUE), n = 5))
    print(head(sort(table(relationship$from), decreasing = TRUE), n = 5))

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
    