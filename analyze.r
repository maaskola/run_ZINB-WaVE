#!/usr/bin/env Rscript

install_package <- FALSE

if (install_package) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("zinbwave")
}

library(zinbwave)
library(scRNAseq)
library(matrixStats)
library(magrittr)
library(ggplot2)
library(biomaRt)

# Register BiocParallel Serial Execution
BiocParallel::register(BiocParallel::SerialParam())

mob_paths <- list.files("~/expression/mob16", "mob_.*.tsv.gz", full.names=TRUE)
brain_paths <- list.files("~/expression/mob-hippocampus/transposed/hippocampus/", "wt_rep.*.tsv.gz", full.names=TRUE)

# assume paths are for count matrices with spots in rows and genes in columns
load_data <- function(paths) {
  sexps <- list()
  for (path in paths) {
    counts <- st.load.matrix(path)
    print(counts[1:10,1:10])
    sexp <- SummarizedExperiment(assays=list(counts=counts))
    # sexp <- SummarizedExperiment(assays=list(counts=t(counts)),
    #                              colData=colData)
                                 # rowRanges=rowRanges, colData=colData)
    sexps[[path]] <- sexp
  }
  sexps
}

# assume paths are for count matrices with spots in rows and genes in columns
load_data2 <- function(paths) {
  counts <- list()
  colData <- c()
  for (path in paths) {
    print(path)
    count <- st.load.matrix(path)
    print(count[1:10,1:10])
    counts[[path]] = count
  }

  genes <- c()
  for (count in counts)
    genes <- union(genes, rownames(count))

  print(str(genes))
  section <- c()
  coord <- c()
  m <- matrix(0, nrow=length(genes), ncol=0)
  rownames(m) <- genes
  idx <- 1
  print(str(m))
  for (count in counts) {
    n <- matrix(0, nrow=length(genes), ncol=ncol(count))
    rownames(n) <- genes
    colnames(n) <- paste(idx, colnames(count))
    n[rownames(count),] = count
    m <- cbind(m, n)
    section <- c(section, rep(idx, ncol(count)))
    coord <- c(coord, colnames(count))
    print(idx)
    print(str(m))
    idx <- idx + 1
  }

  print(str(section))
  print(table(section))
  print(str(coord))

  coord <- matrix(unlist(lapply(strsplit(coord, "x"), as.numeric)), ncol=2, byrow=TRUE)

  colData <- list(section=section, x=coord[,1], y=coord[,2])

  print(m[1:10,1:10])
  sexp <- SummarizedExperiment(assays=list(counts=m), colData=colData)
    # sexp <- SummarizedExperiment(assays=list(counts=t(counts)),
    #                              colData=colData)
                                 # rowRanges=rowRanges, colData=colData)
  sexp
}

join_data <- function(sexps) {
  genes <- c()
  for (sexp in sexps)
    genes <- union(genes, rownames(sexps))
  
}

analyze <- function(expr, K=2, num_genes=100) {
  print("Before filtering")
  print(dim(expr))
  filter <- rowSums(assay(expr)>5)>5
  print("Gene filter statistics")
  print(table(filter))

  expr <- expr[filter,]
  print("After gene filtering")
  print(dim(expr))

  filter <- colSums(assay(expr)>5)>5
  print("Spot filter statistics")
  print(table(filter))

  expr <- expr[,filter]
  print("After spot filtering")
  print(dim(expr))

  assay(expr) %>% log1p %>% rowVars -> vars
  names(vars) <- rownames(expr)
  vars <- sort(vars, decreasing = TRUE)
  # print(head(vars))

  expr <- expr[names(vars)[1:num_genes],]

  assayNames(expr)[1] <- "counts"

  zinb <- zinbwave(expr, K=K, epsilon=1000)

  # print(str(zinb))
  list(expr=expr, zinb=zinb)
}

visualize_it <- function(expr, zinb, output_prefix=NULL) {
  coords <- cbind(x=expr$x, y=expr$y)
  sections <- unique(sort(expr$section))

  W <- reducedDim(zinb)
  # print(str(W))
  if (!is.null(output_prefix)) {
    write.table(W, file=paste0(output_prefix, "ZINB-WaVE-output.tsv"), sep="\t", quote=FALSE)
    write.table(colData(expr), file=paste0(output_prefix, "ZINB-WaVE-colData.tsv"), sep="\t", quote=FALSE)
    for (section_idx in sections) {
      these <- expr$section==section_idx
      w <- W[these,]
      rownames(w) <- gsub(".* ", "", rownames(w))
      write.table(w, file=paste0(output_prefix, "ZINB-WaVE-output-section", section_idx, ".tsv"), sep="\t", quote=FALSE)
    }
  }

  data.frame(W, section=as.character(expr$section)) %>%
    ggplot(aes(W1, W2, colour=section)) + geom_point() + 
    scale_color_brewer(type = "qual", palette = "Set1") +
    theme_classic()
    #scale_color_brewer(type = "qual", palette = "Set1") + theme_classic()
  # data.frame(W, bio=colData(expr)$Biological_Condition,
  #          coverage=colData(expr)$Coverage_Type) %>%
  #   ggplot(aes(W1, W2, colour=bio, shape=coverage)) + geom_point() + 
  #   scale_color_brewer(type = "qual", palette = "Set1") + theme_classic()

  n <- ceiling(sqrt(length(sections)))
  pdf(paste0(output_prefix, "ZINB-WaVE-visual.pdf"), width=n*6, height=n*6)
  # print(str(coords))
  for (col_idx in 1:ncol(W)) {
    par(mfrow=c(n, n))
    for (section_idx in sections) {
      these <- expr$section==section_idx
      w <- W[these,col_idx]
      # print(str(w))
      visualize(w, coords=coords[these,])
    }
  }
  dev.off()
}

doit <- function(expr, K=2, output_prefix="./", num_genes=100) {
  res <- analyze(expr, K=K, num_genes=num_genes)
  visualize_it(res$expr, res$zinb, output_prefix=output_prefix)
  res
}

main_fluidigm <- function() {
  data("fluidigm")
  doit(fluidigm)
}

main <- function() {
  expr <- load_data2(mob_paths[1:4])
  doit(expr)
}

if(!interactive()) {
  main()
}
