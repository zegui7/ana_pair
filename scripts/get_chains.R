library(bio3d)
library(magrittr)

args <- commandArgs(trailingOnly = T)

in_path <- args[1]
cat(in_path)
pdb_files <- list.files(path = in_path,pattern = "*pdb",full.names = T)

get_chains <- function(pdb_file) {
 pdb <- read.pdb(pdb_file)
 pdb$atom$resid[grep(pattern = "HID|HIE|HIP",x = pdb$atom$resid)] <- "HIS"
 pdb$atom$resid[grep(pattern = "CYX",x = pdb$atom$resid)] <- "CYS" 
 pdb$atom$chain <- seq(10,length(pdb$atom$chain) + 9)
 chains <- chain.pdb(pdb,rtn.vec = T,blank = "X")
 write.table(chains,file = "chains")
}

get_chains(pdb_file = pdb_files[1])
