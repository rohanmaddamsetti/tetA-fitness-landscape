## merge-all-growth-curve-data.R by Rohan Maddamsetti.
## This script makes a big dataframe with all the growth curves in DH5a
## or in K-12 with tetA-GFP clones, and writes to file as a CSV file.

library(tidyverse)


## get the DH5a tetA clone growth data.
Jan30.LB.data <- read.csv("../results/growth-results/2022-01-30-LB-Tet0-B30-clones-tidy.csv") %>%
    mutate(Date="2022-01-30")

Jan30.Tet20.data <- read.csv("../results/growth-results/2022-01-30-LB-Tet20-B30-clones-tidy.csv") %>%
    mutate(Date="2022-01-30") 

Jan30.Tet50.data <- read.csv("../results/growth-results/2022-01-30-LB-Tet50-B30-clones-tidy.csv") %>%
    mutate(Date="2022-01-30")

Feb01.LB.data <- read.csv("../results/growth-results/2022-02-01-LB-Tet0-B30-clones-tidy.csv") %>%
    mutate(Date="2022-02-01")

Feb01.Tet20.data <- read.csv("../results/growth-results/2022-02-01-LB-Tet20-B30-clones-tidy.csv") %>%
    mutate(Date="2022-02-01")

Feb01.Tet50.data <- read.csv("../results/growth-results/2022-02-01-LB-Tet50-B30-clones-tidy.csv") %>%
    mutate(Date="2022-02-01")

Feb03.LB.data <- read.csv("../results/growth-results/2022-02-03-LB-Tet0-B30-clones-tidy.csv") %>%
    mutate(Date="2022-02-03")

Feb03.Tet20.data <- read.csv("../results/growth-results/2022-02-03-LB-Tet20-B30-clones-tidy.csv") %>%
    mutate(Date="2022-02-03")

Feb03.Tet50.data <- read.csv("../results/growth-results/2022-02-03-LB-Tet50-B30-clones-tidy.csv") %>%
    mutate(Date="2022-02-03")



## join into a big dataframe for DH5a growth curves.
DH5a.tetA.growth.curve.data <- Jan30.LB.data %>%
    full_join(Jan30.Tet20.data) %>%
    full_join(Jan30.Tet50.data) %>%
    full_join(Feb01.LB.data) %>%
    full_join(Feb01.Tet20.data) %>%
    full_join(Feb01.Tet50.data) %>%
    full_join(Feb03.LB.data) %>%
    full_join(Feb03.Tet20.data) %>%
    full_join(Feb03.Tet50.data) %>%
    ## add shared metadata.
    ## Cultures for this experiment were preconditioned in LB for 24h.
    mutate(Strain="DH5a") %>%
    mutate(preconditioning_length="24h") %>%
    mutate(PreconditioningTet = 0) %>%
    ## make columns consistent with K12 growth curves.
    mutate(Plasmid = fct_recode(Treatment,
                                `No plasmid` = "no_plasmid",
                                p15A = "p15A_plasmid",
                                pUC = "pUC_plasmid")) %>%
    select(-Treatment) %>%
    ## we don't need raw OD600 measurements in the combined dataframe.
    select(-RawOD600) %>%
    ## All of the DH5a pUC clones have pUC, so EvolvedPlasmid == Plasmid.
    ## IMPORTANT TODO: go through my qPCR data and double-check this!
    mutate(lost_pUC = FALSE) %>%
    mutate(EvolvedPlasmid = Plasmid)


## get the barcoded K12 tetA-GFP clone growth data.
K12.tetA.GFP.growth.curve.data <- read.csv("../results/tetA-GFP-growth-results/August2023-tetA-GFP-clone-growth-curves.csv") %>%
    mutate(Strain="K12_barcode") %>%
    ## only one clone sampled from each EvolutionReplicate Population.
    mutate(Clone = 1)


## now merge the all the growth curves,
DH5a.and.K12.tetA.growth.curve.data <- full_join(
    DH5a.growth.curve.data, K12.tetA.GFP.growth.curve.data)
## and write to file.
write.csv(DH5a.and.K12.tetA.growth.curve.data,
          "../results/combined-DH5a-K12-growth-curves.csv", quote=FALSE, row.names=FALSE)
