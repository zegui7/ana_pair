options(warn = -1)
library(bio3d)
library(magrittr)

args <- commandArgs(trailingOnly = T)

pdb_file <- args[1]
pdb_id <- args[2]
out_path <- args[3]

get_fluctuations <- function(pdb_file,pdb_id,out_path) {
 pdb <- read.pdb(pdb_file)
 nma_obj <- nma(pdb = pdb)
 fluctuations <- nma_obj$fluctuations
 out <- data.frame(paste0(pdb$atom$resid[pdb$calpha],
                          pdb$atom$resno[pdb$calpha],
                          pdb$atom$chain[pdb$calpha]),
                   nma_obj$fluctuations)
 write.table(x = out,file = out_path,row.names = F, col.names = F,
             sep = ',',quote = FALSE)
 gc()
}

get_fluctuations(pdb_file,pdb_id, out_path)