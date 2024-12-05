#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Input parameters
params.input_rds = "./data/AD_microglia_data.rds"

// Define the pipeline workflow
workflow {

    input_rds = file(params.input_rds)

    loadData_out = loadData(input_rds)

    de_out = differentialExpression(loadData_out.counts, loadData_out.metadata)

    go_out = geneOntology(de_out.de_results)

    plots = plotResults(loadData_out.counts, de_out.de_results, go_out.go_results)
}

// Process to load Seurat object
process loadData {
    tag "Load Seurat data"
    
    publishDir "./results/loading", mode: 'copy'

    input:
    path input_seurat_rds

    output:
    path "metadata.rds", emit: metadata
    path "counts_matrix.rds", emit: counts

    container = 'kzhang04/de-pipeline:0.1.0'

    script:
    """
    Rscript -e "
    library(Seurat);
    microglia_data <- readRDS('${input_seurat_rds}');
    metadata <- microglia_data@meta.data;
    counts <- GetAssayData(microglia_data, slot = 'counts');
    saveRDS(metadata, 'metadata.rds');
    saveRDS(counts, 'counts_matrix.rds');
    "
    """
}

// Process for differential expression analysis
process differentialExpression {
    tag "Perform DE analysis"

    publishDir "./results/DE", mode: 'copy'

    input:
    path counts
    path metadata

    output:
    path "DE_results.rds", emit: de_results

    container = 'kzhang04/de-pipeline:0.1.0'

    script:
    """
    Rscript -e "
    library(DESeq2);
    counts <- readRDS('$counts');
    metadata <- readRDS('$metadata');
    counts <- as.matrix(counts)+1
    dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ disease);
    dds <- DESeq(dds);
    results <- results(dds);
    saveRDS(results, 'DE_results.rds');
    "
    """
}

// Process for gene ontology analysis
process geneOntology {
    tag "Perform GO analysis"

    publishDir "./results/GO", mode: 'copy'

    input:
    path de_results

    output:
    path "GO_results.rds", emit: go_results

    container = 'kzhang04/de-pipeline:0.1.0'

    script:
    """
    Rscript -e "
    library(clusterProfiler);
    library(org.Hs.eg.db);
    results <- readRDS('$de_results');
    genes <- rownames(results[results\\\$padj < 0.05 & !is.na(results\\\$padj), ]);
    go_results <- enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = 'ENSEMBL', ont = 'BP');
    saveRDS(go_results, 'GO_results.rds');
    "
    """
}

// Process for visualization
process plotResults {
    tag "Plot results"

    publishDir "./results/visualizations", mode: 'copy'

    input:
    path counts
    path de_results
    path go_results

    output:
    path "volcano_plot.png", emit: volcano_plot
    path "heatmap.png", emit: heatmap
    path "go_bar_plot.png", emit: go_plot

    container = 'kzhang04/de-pipeline:0.1.0'

    script:
    """
    Rscript -e "
    library(DESeq2);
    library(ggplot2);
    library(pheatmap);
    library(clusterProfiler);
    library(org.Hs.eg.db);
    
    # Load differential expression results
    de_results <- readRDS('$de_results');
    de_results <- as.data.frame(de_results);
    
    # Volcano Plot
    # Calculate -log10(p-value)
    de_results\\\$log_pval <- -log10(de_results\\\$padj);
    volcano_plot <- ggplot(de_results, aes(x = log2FoldChange, y = log_pval)) +
        geom_point(aes(color = padj < 0.05), alpha = 0.6) +
        scale_color_manual(values = c('black', 'red')) +
        labs(title = 'Volcano Plot of Differentially Expressed Genes', x = 'Log2 Fold Change', y = '-Log10(p-value)') +
        theme_minimal() +
        theme(
            panel.background = element_rect(fill = 'white'),
            plot.background = element_rect(fill = 'white')
        );
    ggsave('volcano_plot.png', volcano_plot);
    
    #Load counts data
    counts <- readRDS('$counts');

    # Heatmap
    # Select top 50 most significant genes (by adjusted p-value)
    top_genes <- de_results[order(de_results\\\$padj), ][1:50, ];
    top_genes_list <- rownames(top_genes);
    heatmap_data <- counts[top_genes_list, , drop = FALSE];
    
    # Normalize counts for heatmap
    heatmap_data <- t(scale(t(heatmap_data)));
    
    # Generate heatmap
    heatmap_plot <- pheatmap(heatmap_data, 
                             cluster_rows = TRUE, 
                             cluster_cols = TRUE, 
                             show_rownames = FALSE, 
                             show_colnames = FALSE, 
                             main = 'Heatmap of Top 50 DE Genes');
    ggsave('heatmap.png', heatmap_plot\\\$gtable);
    
    # GO Enrichment Plot (Bar plot of the top GO terms)
    go_results <- readRDS('$go_results');
    go_results <- as.data.frame(go_results);
    go_plot <- ggplot(go_results[1:10, ], aes(x = reorder(Description, pvalue), y = -log10(pvalue))) +
        geom_bar(stat = 'identity', aes(fill = pvalue)) +
        labs(title = 'Top 10 GO Enrichment Terms', x = 'GO Term', y = '-Log10(p-value)') +
        theme_minimal() +
        theme(
            panel.background = element_rect(fill = 'white'),
            plot.background = element_rect(fill = 'white')
        ) +
        coord_flip();
    ggsave('go_bar_plot.png', go_plot);
    "
    """
}
