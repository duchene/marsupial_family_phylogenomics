library(phangorn)
dats <- grep("phy$", dir(), value = T)
dirs <- list.dirs()
count <- 0
for(i in 1:length(dats)){
      if(count < 100){
      if(!paste0("./", dats[i], ".ML.folder") %in% dirs){
      count <- count + 1
      system(paste0("mkdir ", gsub(".phy", "", dats[i]), ".ML.folder"))
      setwd(paste0(gsub(".phy", "", dats[i]), ".ML.folder"))
      system(paste0("cp ../", dats[i], " ../run.mladeq.Rscript ../run.mladeq.sh ."))
      dat <- read.dna(dats[i])
      rownames(dat) <- gsub(" ", "_", rownames(dat))
      write.dna(dat, file = dats[i])
      system("qsub run.mladeq.sh")
      setwd("..")
      }
      }
}