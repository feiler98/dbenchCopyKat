# R function for running CopyKat

library(copykat)
library(Seurat)

csv_to_seurat <- function(path_rcm) {
  raw_rcm <- read.csv(path_rcm, header = T)
  rownames(raw_rcm) <-  make.names(raw_rcm$Gene, unique = TRUE)
  raw_rcm$Gene <- NULL
  seurat_obj <- CreateSeuratObject(counts = raw_rcm, assay="RNA")
  return (GetAssayData(seurat_obj ,assay = 'RNA', layer = 'counts'))
}

r_run_copykat <- function(path_file,
                          sam.name,
                          n.cores=10,
                          ngene.chr=5,
                          wind.size=25,
                          KS.cut=0.1,
                          LOW.DR=0.05,
                          UP.DR=0.1,
                          norm.cell.names=""
                          ) {

    # get counts from seurat object
    counts <- csv_to_seurat(path_file)

    # pipeline strictly for human data only constructed with hg38 ref genome
    copykat(rawmat=counts,
            id.type = "S",
            ngene.chr=ngene.chr,
            win.size=wind.size,
            KS.cut=KS.cut,
            LOW.DR=LOW.DR,
            UP.DR=UP.DR,
            sam.name=sam.name,
            distance="euclidean",
            norm.cell.names=norm.cell.names,
            output.seg="FALSE",
            plot.genes="TRUE",
            genome="hg20",
            n.cores=n.cores)
}
