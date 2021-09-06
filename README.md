# SeqTracks-plot-R
This is a set of functions for creating figures/plots from seqeuncing data
* `plot_tl`:       Generates reads coverage track from TimeLapse experiment showing 0, 1, 2, 3, 4, 5+ mutation contatining reads.
* `plot_track`:    Gerenerates reads coverage track from single sample or overlay of mnultiple samples.
* `plot_BarGene`:  Allows plotting of descrete scale data onto continuous scale annotation.


## plot_tl
```R
plot_tl (dir, name, chr, start, end, yscale = NULL, labx = TRUE)
```
    dir     -   directory to serch for .bigWig files
    name    -   a unique portion of the file names that selects all mutation track files for both strands e.g. "NMD_rep1"
    chr     -   chromosome identifier for plotting window
    start   -   start (leftmost) position for plotting window <int>
    end     -   end (rightmost) position for plotting window <int>
    yscale  -   optional: vector in format c(min,max) for y axis scale (default: autoscale)
    labx    -   optional: show genomic coordinates on x axis (default: TRUE)
    
Example:

     plot_tl("TracksBigwig/Control/", "1hr_1", "chr6", 36592363, 36607600)
<img src="https://github.com/Mashin6/SeqTracks-plot-R/blob/main/images/plot_tl.png" width="450">



## plot_tracks
```R 
plot_track (dir, name, chr, start, end, yscale = NULL, labx = TRUE, snames = name) 
```
    dir     -   directory to serch for bigWig files
    name    -   either single or vector of names to distinguish files  e.g. "NMD_rep1"
    chr     -   chromosome identifier for plotting window
    start   -   start (leftmost) position for plotting window <int>
    end     -   end (rightmost) position for plotting window <int>
    yscale  -   optional: vector in format c(min,max) for y axis scale (default: autoscale)
    labx    -   optional: show genomic coordinates on x axis (default: TRUE)
    snames  -   optional: sample names for leged (default: file name)

Example:

    plot_track("TracksBigwig/Control", c("1hr_1.TC.0.pos", "1hr_1.TC.1.pos"), "chr6", 36592363, 36607600)
<img src="https://github.com/Mashin6/SeqTracks-plot-R/blob/main/images/plot_track.png" width="450">



## plot_BarGene
```R
plot_BarGene (name_of_gene, annotation, bar_ylab = "", labx = TRUE)
```
    name_of_gene     -   name of gene for plotting, has to match gene_name column values in annotation
    bar_ylab         -   y axis label for barplot (default: off)
    labx             -   show genomic coordinates on x axis (default: TRUE)
    annotation       -   data.frame used for plotting
                         input is dataframe of exon data from ideally 1 but possibly multiple transcripts
        annotation data.frame has to contain these columns: seqnames, start, end, strand, score, sd, transcript_id, gene_name
        each row has information about a single exon: seqnames      - chromosome name
                                                      start         - leftmost edge of exon
                                                      end           - rightmost edge of exon
                                                      strand        - strand information
                                                      score         - value for barplot
                                                      sd            - value for error bars
                                                      transcript_id - transcript identifier this exon belongs to
                                                      gene_nname    - name of gene this exon belongs to

Example:

    plot_BarGene("SRSF3", annot, bar_ylab = "Data +/- SD")
<img src="https://github.com/Mashin6/SeqTracks-plot-R/blob/main/images/plot_BarGene.png" width="350">



## Stacking multiple tracks together
``` R
library(patchwork)
p1 <- plot_BarGene("SRSF3", annot, bar_ylab = "Data +/- SD", labx = FALSE)
p2 <- plot_track("TracksBigwig/Control", c("1hr_1.TC.0.pos", "1hr_1.TC.1.pos"), "chr6", 36593781, 36606397, labx = FALSE, snames = c("Sample1", "Sample2"))
p3 <- plot_tl("TracksBigwig/Control/", "1hr_1", "chr6", 36593781, 36606397) + ylab("Sample")

p1 + p2 + p3 + plot_layout(ncol = 1, heights = c(4, 1, 1),  guides = 'collect')
```
<img src="https://github.com/Mashin6/SeqTracks-plot-R/blob/main/images/plot_combined.png" width="450">




