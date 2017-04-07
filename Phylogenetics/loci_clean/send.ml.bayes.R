library(phangorn)
source("../R/run.mb.R")
load("locinames.Rdata")
birdat <- locinames
dirs <- list.dirs()
count <- 0

for(i in 1:length(birdat)){

if(count < 100){
      if(!paste0("./", gsub("[.]phy", "", birdat[i]), ".Lboot.folder") %in% dirs){
      count <- count + 1

      #dirmlgtr <- paste0(gsub("[.]phy", "", birdat[i]), ".ML.folder")
      #system(paste0("mkdir ", dirmlgtr))
      #setwd(dirmlgtr)
      #system(paste0("cp ../../run.gtr.phyml.sh ../../run.gtr.phyml.Rscript ../", birdat[i], " ."))
      #system("qsub run.gtr.phyml.sh")
      #setwd("..")

      dirmlgtr <- paste0(gsub("[.]phy", "", birdat[i]), ".Lboot.folder")
      system(paste0("mkdir ", dirmlgtr))
      setwd(dirmlgtr)
      system(paste0("cp ../../run.gtr.Lboot.sh ../../run.gtr.Lboot.Rscript ../", birdat[i], " ."))
      system("qsub run.gtr.Lboot.sh")
      setwd("..")

      #dirbayesgtr <- paste0(gsub("[.]phy", "", birdat[i]), ".bayes.folder")
      #system(paste0("mkdir ", dirbayesgtr))
      #setwd(dirbayesgtr)
      #system(paste0("cp ../../run.gtr.mb.sh ../../run.gtr.mb.Rscript ../", birdat[i], " ."))
      #run.mb(birdat[i])
      #system(paste("cp", grep("[.]nex", dir(), value = T), "analysis.nex"))
      #system("qsub run.gtr.mb.sh")
      #setwd("..")
      }
}
}