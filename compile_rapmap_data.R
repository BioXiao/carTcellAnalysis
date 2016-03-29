library(CoverageView)
library(GenomicRanges)
library(rtracklayer)
library(Rsamtools)
library(GenomicAlignments)
library(ggbio)

bam_file <- "~/Box Sync/data/projects/carTcellAnalysis/results/rapmap/150910_D00565_0096_AC6REVANXX/P89-10/lib9207_C6REVANXX.bam"
cov_bam_file <- CoverageView::CoverageBamFile(bam_file)

temp <- readGAlignments(bam_file)
cvg <- coverage(temp)

car_cov <- cvg[["CAR-1"]]
car_smooth_cov <- round(runmean(car_cov, 101, endrule = "constant"))
df <- data_frame(cov = as.vector(window(car_smooth_cov, 1, length(car_cov))))
df %>% ggplot(aes(x = x, y = y)) +
    geom_line() +
    geom_smooth(se = FALSE)

