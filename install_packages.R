packages <- c("Rtools","tidyverse", "cowplot", "ggplot2", "patchwork","dplyr", "fcoex", "ggpubr","SingleR",
              "glmGamPoi", "sctransform", "pheatmap", "Signac", "EnsDb.Mmusculus.v79", "LOLA", "Biostrings", "JASPAR2020", "MotifDb", "BSgenoma.Mmusculus.UCSC.mm10","rtracklayer"
              , "BSgenome", "PMEnrich","PMEnrich", "Mmusculus.background", "glmGamPoi", "DESeq2", "ggrepel")

for (pkg in packages) {
  if(!requireNamespace(pkg, quietly = TRUE)){
    install.packages(pkg, dependencies = TRUE)
  }
}
