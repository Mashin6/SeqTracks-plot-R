#########
# Author: Martin Machyna
# Date: 08/09/21
# Description: Functions for creating sequencing tracks plots and gene annotated barplots in R
#########

plot_tl <- function(dir, name, chr, start, end, yscale = NULL, labx = TRUE) {
    # dir     -   directory to serch for bigWig files
    # name   -   a unique portion of the file names that selects all mutation track files for both strands e.g. "NMD_rep1"
    # chr     -   chromosome identifier for plotting window
    # start   -   start (leftmost) position for plotting window <int>
    # end     -   end (rightmost) position for plotting window <int>
    # yscale  -   optional: vector in format c(min,max) for y axis scale (default: autoscale)
    # labx    -   optional: show genomic coordinates on x axis (default: TRUE)
    #
    # chr, start, end can also be supplied as vectors of the same length
    # Requires: tidyverse, rtracklayer, cowplot, GenomicRanges

    # load packages
    library(tidyverse)
    library(rtracklayer)
    library(cowplot)
    library(GenomicRanges)

    # Get paths to all .bigWig files
    bwFiles <- dir(dir, pattern = paste0(".*", name, ".*bigWig"), full.names = TRUE)
    if(length(bwFiles) == 0) {stop("No files found at provided path")} 


    # Create region that will be plotted
    selectRange <- GRanges(seqnames = chr, ranges = IRanges(start = start, end = end))

    # Load data from .bigWig files
    bwTracks <- NULL
    for (f in bwFiles) {
        n <- basename(f) %>% str_extract('[0-5].(pos|min)')
        bwTracks[[n]] <- import(f, which = selectRange) %>% as.data.frame()
    }

    bwTracks <- bind_rows(bwTracks, .id = "id") %>%                         # combine tracks together and mark them with name
                    mutate(id = str_extract(id, '[0-5]')) %>%               # extract number of mutations from name
                    mutate(start = start - 1) %>%                           # adjust start becuase it is 1-based
                    mutate(score_m = if_else(score < 0, score, 0),          # split minus and plus score
                           score = if_else(score >= 0, score, 0))               
    if(nrow(bwTracks) == 0) {stop("No data retrieved from .bigWig files")} 


    tlPlot <- bwTracks %>% 
                mutate(id = factor(id, levels = c(0:5))) %>%            # assign factors to mutation counts in order to preserve color scale
                ggplot() + 
                    ggplot2::geom_rect(aes(xmin = start, xmax = end, ymin = score_m, ymax = score, fill = id)) +
                    scale_fill_manual(values = c("#c8c8c8", "#fa9696", "#fa0000", "#960000", "#640000", "#320000"),
                                      name = "Muts",
                                      limits = as.character(0:5) ) +
                    geom_hline(yintercept = 0, color = "grey80", size = 0.3)+
                    theme_half_open() +
                    theme(axis.line.x = element_blank(),
                          axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1)  ) 
                    

    if (!labx) {
        tlPlot <- tlPlot +
                    theme(axis.text.x = element_blank(),
                          axis.ticks.x = element_blank())
    }

    if (!is.null(yscale)){
        if (is.na(yscale[1])) {yscale[1] <- ceiling(min(bwTracks$score_m))}
        if (is.na(yscale[2])) {yscale[2] <- floor(max(bwTracks$score))}

        tlPlot <- tlPlot + 
                    coord_cartesian(ylim = c(yscale[1], yscale[2])) +
                    scale_y_continuous(breaks = c(yscale[1], yscale[2]), expand = c(0,0))
    } else {
        tlPlot <- tlPlot +
                    scale_y_continuous(breaks = c(ceiling(min(bwTracks$score_m)), floor(max(bwTracks$score))), expand = c(0,0))
    }

    tlPlot
}




plot_track <- function(dir, name, chr, start, end, yscale = NULL, labx = TRUE, snames = name) {
    # dir     -   directory to serch for bigWig files
    # name    -   either single or vector of names to distinguish files  e.g. "NMD_rep1"
    # chr     -   chromosome identifier for plotting window
    # start   -   start (leftmost) position for plotting window <int>
    # end     -   end (rightmost) position for plotting window <int>
    # yscale  -   optional: vector in format c(min,max) for y axis scale (default: autoscale)
    # labx    -   optional: show genomic coordinates on x axis (default: TRUE)
    # snames  -   optional: sample names for leged (default: file name)
    #
    # chr, start, end can also be supplied as vectors of the same length
    # Requires: tidyverse, rtracklayer, cowplot, GenomicRanges

    # load packages
    library(tidyverse)
    library(rtracklayer)
    library(cowplot)
    library(GenomicRanges)

    # Get paths to all .bigWig files
    bwFiles <- NULL
    for (n in name) {
        bwFiles <- c(bwFiles, dir(dir, pattern = paste0(".*", n, ".*bigWig"), full.names = TRUE))
    }
    if(length(bwFiles) == 0) {stop("No files found at provided path")} 
    

    # Create region that will be plotted
    selectRange <- GRanges(seqnames = chr, ranges = IRanges(start = start, end = end))

    bwTracks <- NULL
    # Load data from .bigWig files
    for (f in bwFiles) {
        n <- name[f == bwFiles]
        bwTracks[[n]] <- import(f, which = selectRange) %>% as.data.frame()
    }

    bwTracks <- bind_rows(bwTracks, .id = "id") %>%                         # combine tracks together and mark them with name
                    mutate(start = start - 1) %>%                           # adjust start becuase it is 1-based
                    mutate(score_m = if_else(score < 0, score, 0),          # split minus and plus score
                           score = if_else(score >= 0, score, 0)) %>%
                    mutate(id = factor(id))                                 # assign factors to file data
    if(nrow(bwTracks) == 0) {stop("No data retrieved from .bigWig files")}               


    tlPlot <- bwTracks %>% 
                ggplot() + 
                    ggplot2::geom_rect(aes(xmin = start, xmax = end, ymin = score_m, ymax = score, fill = id)) +
                    scale_fill_grey(labels = snames,
                                    limits = levels(bwTracks$id)) +
                    geom_hline(yintercept = 0, color = "grey80", size = 0.3)+
                    theme_half_open() +
                    theme(axis.line.x = element_blank(),
                          axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1),
                          legend.title = element_blank())


    if (length(name) == 1) { tlPlot <- tlPlot + theme(legend.position = "none") }

    if (labx) {
        tlPlot <- tlPlot +
                    xlab(bwTracks$seqnames[1])
    }
    else {
        tlPlot <- tlPlot +
                    theme(axis.text.x = element_blank(),
                          axis.ticks.x = element_blank())
    }

    if (!is.null(yscale)){
        if (is.na(yscale[1])) {yscale[1] <- ceiling(min(bwTracks$score_m))}
        if (is.na(yscale[2])) {yscale[2] <- floor(max(bwTracks$score))}

        tlPlot <- tlPlot + 
                    coord_cartesian(ylim = c(yscale[1], yscale[2])) +
                    scale_y_continuous(breaks = c(yscale[1], yscale[2]), expand = c(0,0))
    } else {
        tlPlot <- tlPlot +
                    scale_y_continuous(breaks = c(ceiling(min(bwTracks$score_m)), floor(max(bwTracks$score))), expand = c(0,0))
    }

    tlPlot
}




plot_BarGene <- function(name_of_gene, annotation, bar_ylab = "", labx = TRUE) {
    require(tidyverse)
    require(ggbio)
    require(rtracklayer)
    require(GenomicRanges)
    require(cowplot)
    require(patchwork)

    # name_of_gene     -   name of gene for plotting, has to match gene_name column values in annotation
    # bar_ylab         -   y axis label for barplot (default: off)
    # labx             -   show genomic coordinates on x axis (default: TRUE)
    # annotation       -   data.frame used for plotting
    #                       input is dataframe of exon data from ideally 1 but possibly multiple transcripts
    #     annotation data.frame has to contain these columns: seqnames, start, end, strand, score, sd, transcript_id, gene_name
    #     each row has information about a single exon: seqnames      - chromosome name
    #                                                   start         - leftmost edge of exon
    #                                                   end           - rightmost edge of exon
    #                                                   strand        - strand information
    #                                                   score         - value for barplot
    #                                                   sd            - value for error bars
    #                                                   transcript_id - transcript identifier this exon belongs to
    #                                                   gene_nname    - name of gene this exon belongs to

    # Requires: tidyverse, rtracklayer, ggbio, cowplot, GenomicRanges, patchwork


    annotation_GR <- annotation %>%
                        filter(gene_name == name_of_gene) %>% 
                        makeGRangesFromDataFrame(keep.extra.columns = TRUE)                    

    gene_width <- max(annotation$end) - min(annotation$start)
    bar_width <-  gene_width / nrow(annotation) # range diveded by number of exonic parts

    annotation <- annotation %>% 
                    arrange(start) %>%
                    mutate(bar_start = min(annotation$start) + (bar_width * 1:nrow(annotation)) ) %>% 
                    mutate(score_m = if_else(score < 0, score, 0),          # split minus and plus score
                           score_p = if_else(score >= 0, score, 0))
    

    yscale <- ceiling(max(abs(annotation$score)))

    xlim_min <- min(annotation$start) - gene_width * 0.005
    xlim_max <- max(annotation$end) + gene_width * 0.005

    p1 <- annotation %>%
        ggplot() +
            ggplot2::geom_rect(aes(xmin = bar_start, xmax = bar_start - bar_width * 0.8, ymin = score_m, ymax = score_p)) +
             geom_errorbar(aes(x = bar_start - bar_width * 0.4, y = score, ymax = score + sd, ymin = score - sd), 
                           width = bar_width * 0.2) +
            theme_half_open() +
            background_grid(major = "y", 
                            color.major = "grey80") +
            theme(axis.text.x = element_blank(),
                  axis.title.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.line.x = element_blank()) +
            coord_cartesian(ylim = c(-yscale,  yscale),
                            xlim = c( xlim_min, xlim_max)) +
            labs(title = name_of_gene,
                 subtitle = paste0(annotation$seqnames[1], ": ", min(annotation$start), "-", max(annotation$end)),
                 y = bar_ylab)

    p2 <- data.frame( x = annotation$bar_start - bar_width * 0.4, y = rep(1, nrow(annotation)), exon = 1:nrow(annotation)) %>%
                add_row(x = (annotation$start + annotation$end) / 2, y = rep(0, nrow(annotation)), exon = 1:nrow(annotation)) %>%
                ggplot() +
                geom_line(aes(x = x, y = y, group = exon)) +
                coord_cartesian(xlim = c(xlim_min, xlim_max)) +
                theme_void()

    p3 <- annotation_GR %>% ggplot() +
                    geom_alignment(aes(group = transcript_id), group.selfish = TRUE) +
                    theme_half_open() +
                    background_grid(major = "y", 
                            color.major = "grey80") +
                    theme(axis.line.x = element_blank()) +
                    coord_cartesian(xlim = c(xlim_min, xlim_max))
                    labs(title = name_of_gene)

    if (!labx) {
        p3 <- p3 +
                theme(axis.text.x = element_blank(),
                      axis.ticks.x = element_blank())
    }

    # patchwork lib
    wrap_plots(p1 + p2 + p3@ggplot + 
         plot_layout(ncol = 1, heights = c(2, 1, 1), guides = 'collect'))
        
}



plot_tl("TracksBigwig/Control/", "1hr_1", "chr6", 36592363, 36607600)
plot_track("TracksBigwig/Control", c("1hr_1.TC.0.pos", "1hr_1.TC.1.pos"), "chr6", 36592363, 36607600)
plot_BarGene("SRSF3", annot, bar_ylab = "Data +/- SD")


# Plots stacking with patchwork
p1 <- plot_BarGene("SRSF3", annot, bar_ylab = "Data +/- SD", labx = FALSE)
p2 <- plot_track("TracksBigwig/Control", c("1hr_1.TC.0.pos", "1hr_1.TC.1.pos"), "chr6", 36593781, 36606397, labx = FALSE, snames = c("Sample1", "Sample2"))
p3 <- plot_tl("TracksBigwig/Control/", "1hr_1", "chr6", 36593781, 36606397) + ylab("Sample")


p1 + p2 + p3 + plot_layout(ncol = 1, heights = c(4, 1, 1),  guides = 'collect')



