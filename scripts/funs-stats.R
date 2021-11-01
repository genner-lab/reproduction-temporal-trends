#!/usr/bin/env Rscript

# load libs
suppressPackageStartupMessages({
    library("here")
    library("tidyverse")
    library("lubridate")
    library("magrittr")
    library("glue")
    library("ggthemes")
    library("FactoMineR")
    library("broom")
    library("broom.mixed")
    library("MuMIn")
    library("rfishbase")
    library("glmmTMB")
    library("DHARMa")
    library("emmeans")
    library("performance")
    library("vegan")
    library("xtable")
    #renv::install("cmartin/ggConvexHull")
    library("ggConvexHull")
    library("svglite")
    library("see")
    library("patchwork")
    library("qqplotr")
})

### LOAD FUNCTIONS

# function to set fw species to zero, except for EA where salmonids are recorded
rm_fw <- function(df) {
    if("individualCount" %in% colnames(df) == TRUE) {
        df.rm <- df %>% mutate(individualCount=if_else(order=="Cypriniformes" | species=="Perca fluviatilis" | species=="Cottus gobio" | species=="Esox lucius" | species=="Lampetra fluviatilis", 0, individualCount),
                            weightInGrams=if_else(order=="Cypriniformes" | species=="Perca fluviatilis" | species=="Cottus gobio" | species=="Esox lucius" | species=="Lampetra fluviatilis", 0, weightInGrams)) %>%
                        mutate(individualCount=if_else(partnerID!="EA" & order=="Salmoniformes", 0, individualCount), 
                            weightInGrams=if_else(partnerID!="EA" & order=="Salmoniformes", 0, weightInGrams))
    } else if ("nreads" %in% colnames(df) == TRUE) {
        df.rm <- df %>% mutate(nreads=if_else(order=="Cypriniformes" | species=="Perca fluviatilis" | species=="Cottus gobio" | species=="Esox lucius" | species=="Lampetra fluviatilis", 0, nreads)) %>%
                        mutate(nreads=if_else(partnerID!="EA" & order=="Salmoniformes", 0, nreads)) 
    } else stop(writeLines("Column names must have 'nreads' or 'individualCount'."))
return(df.rm)
}


# fun to edit the species for adult analyses
adult_edit <- function(df,collapse) {
    if(collapse==TRUE) {
    df.ed <- df %>% mutate(species=if_else(genus=="Melanogrammus" | genus=="Merlangius","Merlangius/Melanogrammus",species),
        species=if_else(family=="Ammodytidae",family,species),
        species=if_else(species=="Chelon labrosus" | species=="Chelon ramada","Chelon labrosus/ramada",species),
        species=if_else(family=="Triglidae",family,species),
        species=if_else(family=="Triglidae",family,species),
        species=if_else(genus=="Pomatoschistus",genus,species),
        species=if_else(genus=="Gobius",genus,species),
        species=if_else(genus=="Aphia" | genus=="Crystallogobius","Aphia/Crystallogobius",species))
    } else if(collapse==FALSE) {
        df.ed <- df
    } else stop(writeLines("Collapse must be 'TRUE' or 'FALSE'."))
return(df.ed)
}


# fun to edit the species for larvae analyses
larvae_edit <- function(df,collapse) {
    if(collapse==TRUE) {
    df.ed <- df %>% mutate(species=if_else(family=="Lotidae",family,species),
        species=if_else(family=="Ammodytidae",family,species),
        species=if_else(family=="Syngnathidae",family,species),
        species=if_else(family=="Gobiidae",family,species),#species=if_else(genus=="Lebetus",genus,species), # dropped Lebetus to make things more simple
        species=if_else(family=="Clupeidae",family,species),
        species=if_else(family=="Blenniidae",family,species))
    } else if(collapse==FALSE) {
        df.ed <- df
    } else stop(writeLines("Collapse must be 'TRUE' or 'FALSE'."))
return(df.ed)
}


# fun to edit the species for egg analyses
egg_edit <- function(df,collapse) {
    if(collapse==TRUE) {
        df.ed <- df %>% mutate(species=if_else(species!="Sardina pilchardus","Actinopterygii",species))
    } else if(collapse==FALSE) {
        df.ed <- df
    } else stop(writeLines("Collapse must be 'TRUE' or 'FALSE'."))
return(df.ed)
}


# annotate taxonomy
annotate_taxonomy <- function(df) {
    df.with.tax <- df %>% mutate(class=if_else(taxonRank=="species",pull(tax.table,class)[match(species,pull(tax.table,sciNameValid))],"NA"),
        order=if_else(taxonRank=="species",pull(tax.table,order)[match(species,pull(tax.table,sciNameValid))],"NA"),
        family=if_else(taxonRank=="species",pull(tax.table,family)[match(species,pull(tax.table,sciNameValid))],"NA"),
        genus=if_else(taxonRank=="species",pull(tax.table,genus)[match(species,pull(tax.table,sciNameValid))],"NA"),
        class=if_else(taxonRank=="class",pull(tax.table,class)[match(species,pull(tax.table,class))],class),
        class=if_else(taxonRank=="family",pull(tax.table,class)[match(species,pull(tax.table,family))],class),
        order=if_else(taxonRank=="family",pull(tax.table,order)[match(species,pull(tax.table,family))],order),
        family=if_else(taxonRank=="family",pull(tax.table,family)[match(species,pull(tax.table,family))],family),
        class=if_else(taxonRank=="genus",pull(tax.table,class)[match(species,pull(tax.table,genus))],class),
        order=if_else(taxonRank=="genus",pull(tax.table,order)[match(species,pull(tax.table,genus))],order),
        family=if_else(taxonRank=="genus",pull(tax.table,family)[match(species,pull(tax.table,genus))],family),
        genus=if_else(taxonRank=="genus",pull(tax.table,genus)[match(species,pull(tax.table,genus))],genus))
    return(df.with.tax)
}


# function to remove fw species and collapse where no species resolution in either eDNA or trad
# collapse taxonomy 
# rmfw sets the freshwater fish reads to zero (does not remove them)
# lifestage must be "Adult", "Larvae", or "Eggs"
# collapse==FALSE lets us remove the egg/larvae data but keep the names unchanged
collapse_taxonomy <- function(df,rmfw,lifestage,collapse){
    # annotate taxonomy
    df.with.tax <- df %>% annotate_taxonomy()
# remove the freshwater
    if(rmfw==TRUE) {
        df.with.tax %<>% rm_fw()
    } else if(rmfw==FALSE) {
        df.with.tax <- df.with.tax
    } else stop(writeLines("rmfw must be 'TRUE' or 'FALSE'."))
# to remove unwanted lifestages from trad data
    if("individualCount" %in% colnames(df.with.tax) == TRUE) {
        df.with.tax %<>% filter(lifeStage==lifestage) %>% rename(lifestage=lifeStage)
    } else if ("nreads" %in% colnames(df.with.tax) == TRUE) {
        df.with.tax %<>% mutate(lifestage=lifestage)
    } else stop(writeLines("Column names must have 'nreads' or 'individualCount'."))
# to collapse the names
    if(lifestage=="Adult") {
        df.with.tax %<>% adult_edit(collapse=collapse)
    } else if (lifestage=="Eggs") {
        df.with.tax %<>% filter(class=="Actinopterygii") %>% egg_edit(collapse=collapse)
    } else if (lifestage=="Larvae") {
        df.with.tax %<>% adult_edit(collapse=collapse) %>% larvae_edit(collapse=collapse)
#    } else if (lifestage=="Traits") {
#        df.with.tax %<>% traits_edit()
    } else stop(writeLines("Lifestage must be 'Adult', 'Eggs', or 'Larvae'."))
    df.with.tax %<>% select(-taxonRank,-class,-order,-family,-genus)
return(df.with.tax)
}


# function to create grouping variables
add_groups_vars <- function(df,spatialgroup,temporalgroup) {
    completed.df <- df
        if(temporalgroup=="month") {
            completed.df %<>% mutate(temporalGroup=paste(year,month,"1",sep="-")) } 
        else if(temporalgroup=="day"){
            completed.df %<>% mutate(temporalGroup=paste(year,month,day,sep="-")) } 
        else stop(writeLines("temporalgroup must be 'month' or 'day'."))
        if(spatialgroup=="locality") {
            completed.df %<>% mutate(spatialGroup=localityID) } 
        else if(spatialgroup=="partner") {
            completed.df %<>% mutate(spatialGroup=partnerID) } 
        else if(spatialgroup=="site"){
            completed.df %<>% mutate(spatialGroup=localitySite) } 
        else stop(writeLines("spatialgroup must be 'locality', 'site', or 'partner'."))
    return(completed.df)
}


# function to filter by number of no-catch events
catch_filter <- function(df,filter,maxprop0){
        df.f <- df %>% filter(!is.na(count.edna) & !is.na(count.traditional))
        if(filter=="either0") {
            df.new <- df.f %>% group_by(species) %>% mutate(zeroEvent=if_else(count.edna==0 | count.traditional==0,TRUE,FALSE),nEvents=n(),nZeroEvents=length(which(zeroEvent==TRUE)),propZeroEvents=nZeroEvents/nEvents) %>% ungroup()
        } else if(filter=="both0") {
            df.new <- df.f %>% group_by(species) %>% mutate(zeroEvent=if_else(count.edna==0 & count.traditional==0,TRUE,FALSE),nEvents=n(),nZeroEvents=length(which(zeroEvent==TRUE)),propZeroEvents=nZeroEvents/nEvents) %>% ungroup()
        } else stop(writeLines("filter must be: 'either0' or 'both0'."))
        df.new %<>% filter(propZeroEvents <= maxprop0) %>% 
        select(-zeroEvent,-nEvents,-nZeroEvents)
    return(df.new)
}


# log plus1 fun
logplus1 <- function(x){
    xt <- log(x+1)
    return(xt)
}


# 4th root function
fourth_root <- function(x){
    xt <- x^(1/4)
    return(xt)
}


# standard error function
se <- function(x) sqrt(var(x)/length(x))


# fun to replace NaNs with NAs
# see also mutate_all(~replace(.,is.nan(.),NA)) 
replace_nan <- function(x){
    y <- ifelse(is.nan(x),NA,x)
    return(y)
}


# fun to get life history data from fishbase
annotate_fishbase <- function(df,collapse) {
    # names
    spp.names <- df %>% pull(sciNameValid)
    # spawning
    spawn.tab <- spawning(spp.names,server="fishbase") %>% 
        select(Species,Spawningarea,Jan,Feb,Mar,Apr,May,Jun,Jul,Aug,Sep,Oct,Nov,Dec) %>% 
        filter(!is.na(Species)) %>% 
        rename(sciNameValid=Species,spawningArea=Spawningarea) %>% 
        pivot_longer(cols=c(-sciNameValid,-spawningArea),names_to="month",values_to="spawning") %>% 
        mutate(spawning=if_else(spawning==111,TRUE,NA,FALSE)) %>% 
        group_by(sciNameValid) %>% 
        mutate(hasSpawnData=if_else(TRUE %in% spawning,TRUE,FALSE)) %>%
        ungroup() %>% 
        group_by(sciNameValid,hasSpawnData,month) %>%
        mutate(nSpawn=sum(spawning==TRUE,na.rm=TRUE),nTot=length(month),propSpawningByMonth=nSpawn/nTot) %>%
        ungroup() %>% 
        distinct(sciNameValid,hasSpawnData,month,propSpawningByMonth) %>% 
        mutate(propSpawningByMonth=if_else(hasSpawnData==FALSE,NaN,propSpawningByMonth)) %>% 
        mutate(spawningByMonth=if_else(propSpawningByMonth==0,FALSE,TRUE,NA)) %>% 
        select(-hasSpawnData)
    # custom spawning
    custom.tab <- read_csv(file=here("assets/fishbase-spawning-custom.csv"),show_col_types=FALSE) %>% 
        dplyr::select(-notes) %>%
        filter(!is.na(lifestyle)) %>% 
        rename(sciNameValid=species) %>%
        pivot_longer(cols=c(-sciNameValid,-lifestyle),names_to="month",values_to="spawning") %>% 
        mutate(spawning=if_else(spawning==111,TRUE,NA,FALSE)) %>% 
        rename(spawningByMonthUK=spawning)
    # larval plankton period
    larvae.tab <- larvae(spp.names,server="fishbase") %>% 
        filter(!is.na(Species)) %>% 
        select(Species, ends_with("Larv")) %>% 
        group_by(Species) %>% 
        slice_tail(n=1) %>% 
        rename(sciNameValid=Species) %>% 
        rename_with(~str_replace(.,"Larv","")) %>% 
        pivot_longer(cols=c(-sciNameValid),names_to="month",values_to="larvaeByMonth") %>% #filter(sciNameValid=="Perca fluviatilis")
        mutate(larvaeByMonth=if_else(larvaeByMonth==-1,TRUE,FALSE,NA)) %>%
        group_by(sciNameValid) %>% 
        mutate(hasLarvaeData=if_else(TRUE %in% larvaeByMonth,TRUE,FALSE)) %>%
        ungroup() %>% 
        mutate(larvaeByMonth=if_else(hasLarvaeData==FALSE,NA,larvaeByMonth)) %>% 
        select(-hasLarvaeData)
    # fecundity
    fec.tab <- fecundity(spp.names,server="fishbase") %>% 
        select(Species,Locality,FecundityMin,FecundityMax) %>% 
        filter(!is.na(Species)) %>% 
        rename(sciNameValid=Species,locality=Locality,fecundityMin=FecundityMin,fecundityMax=FecundityMax) %>%
        group_by(sciNameValid) %>%
        mutate(fecundityMinMean=suppressWarnings(mean(fecundityMin,na.rm=TRUE)),fecundityMaxMean=suppressWarnings(mean(fecundityMax,na.rm=TRUE)),fecundityMean=suppressWarnings(mean(c(fecundityMinMean,fecundityMaxMean),na.rm=TRUE))) %>%
        ungroup() %>% 
        distinct(sciNameValid,fecundityMean)
    # lifestyle/species
    life.tab <- rfishbase::species(spp.names,server="fishbase") %>% 
        filter(!is.na(Species)) %>% 
        select(Species,DemersPelag,Length) %>% 
        #mutate(DemersPelag=str_replace_all(DemersPelag,"bathydemersal","demersal")) %>%
        rename(sciNameValid=Species,lifestyleFishbase=DemersPelag,length=Length)
    # length
    length.tab <- length_weight(spp.names,server="fishbase") %>% 
        filter(!is.na(Species)) %>% 
        select(Species,LengthMax,Type) %>% 
        rename(sciNameValid=Species) %>%
        group_by(sciNameValid) %>% 
        summarise(meanMaxLength=mean(LengthMax,na.rm=TRUE),.groups="drop")
    # ecology
    ecology.tab <- ecology(spp.names,server="fishbase") %>% 
        filter(!is.na(Species)) %>% 
        select(Species,DietTroph,FoodTroph,DietTLu,Shoaling,Schooling) %>%
        rename(sciNameValid=Species) %>% 
        group_by(sciNameValid) %>% 
        slice_head(n=1) %>% 
        mutate(meanTrophicLevel=mean(c(DietTroph,FoodTroph,DietTLu),na.rm=TRUE), shoaling=if_else(Shoaling==-1 | Schooling==-1,TRUE,FALSE,NA)) %>%
        ungroup() %>%
        select(sciNameValid,meanTrophicLevel,shoaling)
    # join and clean
    df.joined <- df %>% 
        left_join(spawn.tab,by="sciNameValid") %>% 
        left_join(custom.tab,by=c("sciNameValid","month")) %>% 
        left_join(larvae.tab,by=c("sciNameValid","month")) %>% 
        left_join(fec.tab,by="sciNameValid") %>% 
        left_join(life.tab,by="sciNameValid") %>% 
        left_join(length.tab,by="sciNameValid") %>% 
        left_join(ecology.tab,by="sciNameValid")
    df.joined %<>% relocate(lifestyle,.after=lifestyleFishbase) %>%
        mutate_if(is.numeric,list(replace_nan)) %>%
        rename(species=sciNameValid)
#df.joined %<>% mutate(dblBinBreed=if_else(!!as.name(spawndata)==0,0,1),fctBinBreed=factor(if_else(!!as.name(spawndata)==0,"non-breeding","breeding"),levels=c("non-breeding","breeding")))
#mutate(dblBinLarvae=if_else(lglLarvae==TRUE,1,0),fctBinLarvae=factor(if_else(lglLarvae==TRUE,"larvae-present","larvae-not-present"),levels=c("larvae-present","larvae-not-present"))) %>%
    # collapse
    if(collapse=="none") {
        df.joined <- df.joined
    } else if(collapse=="Adult") {
        df.joined %<>%
            filter(!is.na(lifestyle)) %>%
            adult_edit(collapse=TRUE) %>%
            collapse_fishbase()
    } else if(collapse=="Larvae") {
        df.joined %<>%
            filter(!is.na(lifestyle)) %>%
            adult_edit(collapse=TRUE) %>% larvae_edit(collapse=TRUE) %>%
            collapse_fishbase()
    } else if(collapse=="Eggs") {
        df.joined %<>%
            filter(!is.na(lifestyle)) %>%
            egg_edit(collapse=TRUE) %>%
            collapse_fishbase()
    } else stop(writeLines("collapse must be 'none', 'Adult', 'Larvae', or 'Eggs'."))
return(df.joined)
}

# helper to collapse fishbase
collapse_fishbase <- function(df) {
    df.new <- df %>%
            group_by(species) %>%
            mutate(nd=n_distinct(lifestyleFishbase),lifestyleFishbase=if_else(nd>1,lifestyle,lifestyleFishbase)) %>%
            ungroup() %>%
            group_by(class,order,family,species,lifestyleFishbase,lifestyle,month) %>%
            summarise(propSpawningByMonth=max(propSpawningByMonth,na.rm=FALSE),
                spawningByMonth=if_else(TRUE %in% spawningByMonth,TRUE,FALSE),
                spawningByMonthUK=if_else(TRUE %in% spawningByMonthUK,TRUE,FALSE),
                larvaeByMonth=if_else(TRUE %in% larvaeByMonth,TRUE,FALSE),
                fecundityMean=max(fecundityMean,na.rm=FALSE),
                length=max(length,na.rm=FALSE),
                meanMaxLength=max(meanMaxLength,na.rm=FALSE),
                meanTrophicLevel=max(meanTrophicLevel,na.rm=FALSE),
                shoaling=if_else(TRUE %in% shoaling,TRUE,FALSE),
                .groups="drop")
    return(df.new)
}


# function to filter by number of no-catch events
# comparison can be "individuals" or "weights"
# filter can be "both0" (edna AND trad == 0) or "either0" (edna OR trad == 0)
# maxprop is a filter to remove rows above a threshold
zero_catch_filter <- function(df,comparison,filter,print,maxprop0){
    df.f <- df
        if("nReads" %in% colnames(df.f) == TRUE & comparison=="individuals" & filter=="either0") {
            df.new <- df.f %>% group_by(partnerID,species) %>% mutate(zeroEvent=if_else(nReads==0 | individualsByGroup==0,TRUE,FALSE),nEvents=n(),nZeroEvents=length(which(zeroEvent==TRUE)),propZeroEvents=nZeroEvents/nEvents) %>% ungroup()
        } else if("nReads" %in% colnames(df.f) == TRUE & comparison=="individuals" & filter=="both0") {
            df.new <- df.f %>% group_by(partnerID,species) %>% mutate(zeroEvent=if_else(nReads==0 & individualsByGroup==0,TRUE,FALSE),nEvents=n(),nZeroEvents=length(which(zeroEvent==TRUE)),propZeroEvents=nZeroEvents/nEvents) %>% ungroup()
        } else if("nReads" %in% colnames(df.f) == TRUE & comparison=="weights" & filter=="either0") {
            df.new <- df.f %>% group_by(partnerID,species) %>% mutate(zeroEvent=if_else(nReads==0 | weightInGramsByGroup==0,TRUE,FALSE),nEvents=n(),nZeroEvents=length(which(zeroEvent==TRUE)),propZeroEvents=nZeroEvents/nEvents) %>% ungroup()
        } else if("nReads" %in% colnames(df.f) == TRUE & comparison=="weights" & filter=="both0") {
            df.new <- df.f %>% group_by(partnerID,species) %>% mutate(zeroEvent=if_else(nReads==0 & weightInGramsByGroup==0,TRUE,FALSE),nEvents=n(),nZeroEvents=length(which(zeroEvent==TRUE)),propZeroEvents=nZeroEvents/nEvents) %>% ungroup()
        } else if("nReads" %in% colnames(df.f) == FALSE & comparison=="individuals" & filter=="either0") {
            df.new <- df.f %>% group_by(partnerID,species) %>% mutate(zeroEvent=if_else(readsByGroup==0 | individualsByGroup==0,TRUE,FALSE),nEvents=n(),nZeroEvents=length(which(zeroEvent==TRUE)),propZeroEvents=nZeroEvents/nEvents) %>% ungroup()
        } else if("nReads" %in% colnames(df.f) == FALSE & comparison=="individuals" & filter=="both0") {
            df.new <- df.f %>% group_by(partnerID,species) %>% mutate(zeroEvent=if_else(readsByGroup==0 & individualsByGroup==0,TRUE,FALSE),nEvents=n(),nZeroEvents=length(which(zeroEvent==TRUE)),propZeroEvents=nZeroEvents/nEvents) %>% ungroup()
        } else if("nReads" %in% colnames(df.f) == FALSE & comparison=="weights" & filter=="either0") {
            df.new <- df.f %>% group_by(partnerID,species) %>% mutate(zeroEvent=if_else(readsByGroup==0 | weightInGramsByGroup==0,TRUE,FALSE),nEvents=n(),nZeroEvents=length(which(zeroEvent==TRUE)),propZeroEvents=nZeroEvents/nEvents) %>% ungroup()
        } else if("nReads" %in% colnames(df.f) == FALSE & comparison=="weights" & filter=="both0") {
            df.new <- df.f %>% group_by(partnerID,species) %>% mutate(zeroEvent=if_else(readsByGroup==0 & weightInGramsByGroup==0,TRUE,FALSE),nEvents=n(),nZeroEvents=length(which(zeroEvent==TRUE)),propZeroEvents=nZeroEvents/nEvents) %>% ungroup()
        } else stop(writeLines("filter must be 'either0' or 'both0'; comparison must be 'individuals' or 'weights'."))
    df.new %<>% filter(propZeroEvents <= maxprop0) %>% 
        select(-zeroEvent,-nEvents,-nZeroEvents)
    if(print==TRUE) {
        df.new %>% distinct(species,propZeroEvents) %>% arrange(propZeroEvents,species) %>% print(n=Inf)
    } else {}
    return(df.new)
}


# function to rank species and subset the dataframe by date
# rankmethod must be: 'reads', 'individuals', or 'weights'
# n filters by ranks 1:n
# print prints out the full rankings
rank_and_subset <- function(df,rankmethod,print,n){
    # ranks
    df.rank <- df %>% group_by(partnerID,species)
        if(rankmethod=="reads") {
            df.rank %<>% mutate(sppProps=sum(readsByGroup,na.rm=TRUE))
        } else if(rankmethod=="individuals") {
            df.rank %<>% mutate(sppProps=sum(individualsByGroup,na.rm=TRUE))
        } else if(rankmethod=="weights") {
            df.rank %<>% mutate(sppProps=sum(weightInGramsByGroup,na.rm=TRUE))
        } else stop(writeLines("rankmethod must be: 'reads', 'individuals', or 'weights'."))
    df.rank %<>% ungroup() %>% mutate(rank=dense_rank(desc(sppProps)))
    if(print==TRUE) {
        df.rank %>% distinct(partnerID,species,sppProps,rank) %>% arrange(partnerID,rank) %>% print(n=Inf)
    } else {}
        df.rank %<>% select(-sppProps)
        if(n=="all"){
            df.rank <- df.rank
        } else if(is.numeric(n)) {
            df.rank %<>% filter(rank %in% 1:n)
        } else stop(writeLines("n must be: 'all', or a number."))
        #
#        if(dropspatialgroup=="none" & droptemporalgroup=="none"){
#            df.rank <- df.rank
#        } else if(dropspatialgroup=="site" & droptemporalgroup=="day"){
#            df.rank %<>% select(-eventID,-sampleHash,-replicateFilter,-replicatePCR,-nReads,-readsBySampleTotal) %>% distinct()
#        } else if(dropspatialgroup=="site" & droptemporalgroup=="month"){
#            df.rank %<>% select(-eventID,-sampleHash,-replicateFilter,-replicatePCR,-nReads,-readsBySampleTotal,-eventDate,-day) %>% distinct()
#        } else if(dropspatialgroup=="locality" & droptemporalgroup=="day"){
#            df.rank %<>% select(-eventID,-sampleHash,-replicateFilter,-replicatePCR,-nReads,-readsBySampleTotal,-localitySite) %>% distinct()
#        } else if(dropspatialgroup=="locality" & droptemporalgroup=="month"){
#            df.rank %<>% select(-eventID,-sampleHash,-replicateFilter,-replicatePCR,-nReads,-readsBySampleTotal,-localitySite,-eventDate,-day) %>% distinct()
#        } else if(dropspatialgroup=="partner" & droptemporalgroup=="day"){
#            df.rank %<>% select(-eventID,-sampleHash,-replicateFilter,-replicatePCR,-nReads,-readsBySampleTotal,-localitySite,-localityID) %>% distinct()
#        } else if(dropspatialgroup=="partner" & droptemporalgroup=="month"){
#            df.rank %<>% select(-eventID,-sampleHash,-replicateFilter,-replicatePCR,-nReads,-readsBySampleTotal,-localitySite,-localityID,-eventDate,-day) %>% distinct()
#        } else stop(writeLines("dropgroups must be of: none, or 'day/month' (temporal), and 'site/locality/partner' (spatial)."))
    return(df.rank)
}



# function to join the eDNA and trad survey dataframes
join_and_clean <- function(trad,edna){
    trad.edna.comb <- suppressMessages(left_join(edna,trad)) %>% 
        select(partnerID,primerSet,library,eventID,sampleHash,replicateFilter,replicatePCR,eventDate,year,month,day,localityID,localitySite,
                temporalGroup,spatialGroup,species,lifestage,nReads,maxEfficiency,readsBySampleTotal,readsBySampleProportion,readsByGroupTotal,readsByGroup,readsByGroupProportion,
                individualsByGroup,individualsByGroupRate,individualsByGroupProportion,weightInGramsByGroup,weightInGramsByGroupRate,weightInGramsByGroupProportion) %>%
        arrange(partnerID,primerSet,library,eventID,replicateFilter,replicatePCR) %>% 
        replace_na(list(individualsByGroup=0,individualsByGroupRate=0,individualsByGroupProportion=0,weightInGramsByGroup=0,weightInGramsByGroupRate=0,weightInGramsByGroupProportion=0)) %>%
        mutate(weightInGramsByGroup=if_else(partnerID=="EA" | lifestage!="Adult",NaN,weightInGramsByGroup),
            weightInGramsByGroupRate=if_else(partnerID=="EA" | lifestage!="Adult",NaN,weightInGramsByGroupRate),
            weightInGramsByGroupProportion=if_else(partnerID=="EA" | lifestage!="Adult",NaN,weightInGramsByGroupProportion),
            individualsByGroup=if_else(partnerID!="MBA" & lifestage!="Adult",NaN,individualsByGroup),
            individualsByGroupRate=if_else(partnerID!="MBA" & lifestage!="Adult",NaN,individualsByGroupRate),
            individualsByGroupProportion=if_else(partnerID!="MBA" & lifestage!="Adult",NaN,individualsByGroupProportion),
            individualsByGroup=if_else(localityID=="WHIT" & lifestage!="Adult",NaN,individualsByGroup),
            individualsByGroupRate=if_else(localityID=="WHIT" & lifestage!="Adult",NaN,individualsByGroupRate),
            individualsByGroupProportion=if_else(localityID=="WHIT" & lifestage!="Adult",NaN,individualsByGroupProportion)) %>%
        mutate_all(~replace(.,is.nan(.),NA))
    return(trad.edna.comb)
}


# function to process and collapse the trad and eDNA data
# expand and summarise survey data.
# sppdf is dataframe with merged species.
# eventdf is events dataframe.
# method must be "traditional" or "edna"
# spatialgroup must be "site" (e.g. "L4-bottom") or "locality" (e.g. "L4") or "partner" (e.g. "MBA").
# temporalgroup must be "day" or "month"
# primerdf is the primer bias dataframe
# correct (TRUE/FALSE) to correct the reads by the bias correction factor
expand_and_summarise <- function(df,sppdf,eventdf,primerdf,method,spatialgroup,temporalgroup,correct){
    events.sub <- eventdf %>% distinct(eventID,localityID,localitySite,year,month,day)
    # trad
    if(method=="traditional") {
        completed.df <- df %>% 
            group_by(lifestage,partnerID,eventID,eventDate,fieldNumber,species) %>% 
            summarise(individualCount=sum(individualCount),weightInGrams=sum(weightInGrams),.groups="drop") %>% # merge reads back by species after changing spp names with collapse_taxonomy()
            group_by(partnerID) %>%
            tidyr::complete(tidyr::nesting(lifestage,eventID,eventDate,fieldNumber),species,fill=list(individualCount=0,weightInGrams=0)) %>% # expand with zeros for all species by partnerID. Nesting allows only those combinations of samples present in the dataset.
            ungroup() %>%
            mutate(localityID=pull(events.sub,localityID)[match(eventID,pull(events.sub,eventID))],# annotate with location date etc
                localitySite=pull(events.sub,localitySite)[match(eventID,pull(events.sub,eventID))],
                year=pull(events.sub,year)[match(eventID,pull(events.sub,eventID))],
                month=pull(events.sub,month)[match(eventID,pull(events.sub,eventID))],
                day=pull(events.sub,day)[match(eventID,pull(events.sub,eventID))]) %>% 
            add_groups_vars(spatialgroup=spatialgroup,temporalgroup=temporalgroup) %>% # create grouping variables
            group_by(partnerID,spatialGroup,temporalGroup) %>%  
            mutate(individualsByGroupTotal=sum(individualCount,na.rm=TRUE),weightInGramsByGroupTotal=sum(weightInGrams)) %>% # get event total of counts/weights by grouping variables
            ungroup() %>% 
            group_by(lifestage,partnerID,spatialGroup,temporalGroup,individualsByGroupTotal,weightInGramsByGroupTotal,species) %>%
            summarise(individualsByGroup=sum(individualCount,na.rm=TRUE), # get species totals/rates of counts/weights by grouping variables 
                individualsByGroupRate=mean(individualCount,na.rm=TRUE),
                weightInGramsByGroup=sum(weightInGrams),
                weightInGramsByGroupRate=mean(weightInGrams,na.rm=TRUE),
                .groups="drop") %>%
            mutate(individualsByGroupProportion=individualsByGroup/individualsByGroupTotal,
                weightInGramsByGroupProportion=weightInGramsByGroup/weightInGramsByGroupTotal) %>% # get species proportions of counts/weights by grouping variables 
            mutate(weightInGramsByGroup=if_else(partnerID=="EA",NaN,weightInGramsByGroup),# convert all EA weights to NAs
                weightInGramsByGroupTotal=if_else(partnerID=="EA",NaN,weightInGramsByGroupTotal),
                weightInGramsByGroupRate=if_else(partnerID=="EA",NaN,weightInGramsByGroupRate),
                weightInGramsByGroupProportion=if_else(partnerID=="EA",NaN,weightInGramsByGroupProportion)) %>% 
            mutate_all(~replace(.,is.nan(.),NA))# clean NaNs
    # eDNA
    } else if(method=="edna") {
        completed.df <- df %>% 
            suppressMessages(full_join(sppdf)) %>% # add the names from trad surveys ,by="species"
            group_by(lifestage,primerSet,library,partnerID,eventID,sampleHash,replicateFilter,replicatePCR,eventDate,species) %>% 
            summarise(nreads=sum(nreads,na.rm=TRUE),.groups="drop") %>% # merge reads back by species after changing spp names with collapse_taxonomy()
            group_by(partnerID) %>%
            tidyr::complete(tidyr::nesting(lifestage,primerSet,library,eventID,sampleHash,replicateFilter,replicatePCR,eventDate),species,fill=list(nreads=0)) %>% # expand with zeros for all species by partnerID. Nesting allows only those combinations of samples present in the dataset.
            ungroup() %>%
            add_primer_bias(primerdf=primer.bias,correct=correct) %>%
            filter(!is.na(eventDate)) %>% # remove blanks from the trad species
            mutate(localityID=pull(events.sub,localityID)[match(eventID,pull(events.sub,eventID))], # annotate with location date etc
                localitySite=pull(events.sub,localitySite)[match(eventID,pull(events.sub,eventID))],
                year=pull(events.sub,year)[match(eventID,pull(events.sub,eventID))],
                month=pull(events.sub,month)[match(eventID,pull(events.sub,eventID))],
                day=pull(events.sub,day)[match(eventID,pull(events.sub,eventID))]) %>%
            add_groups_vars(spatialgroup=spatialgroup,temporalgroup=temporalgroup) %>% # create grouping variables
            rename(nReads=nreads) %>%
            group_by(sampleHash) %>%
            mutate(readsBySampleTotal=sum(nReads,na.rm=TRUE)) %>% # get sample total of reads
            ungroup() %>% 
            mutate(readsBySampleProportion=nReads/readsBySampleTotal) %>% #get sample prop of reads
            group_by(partnerID,temporalGroup,spatialGroup) %>%
            mutate(readsByGroupTotal=sum(nReads,na.rm=TRUE)) %>% # get event total of reads by grouping variables
            group_by(partnerID,temporalGroup,spatialGroup,species) %>%
            mutate(readsByGroup=sum(nReads,na.rm=TRUE)) %>% # get species total of reads by grouping variables 
            ungroup() %>%
            mutate(readsByGroupProportion=readsByGroup/readsByGroupTotal) # get species proportion of reads by grouping variables (per species)
    } else stop(writeLines("Method must be 'traditional' or 'edna'."))
    return(completed.df)
}


# add primer bias
add_primer_bias <- function(df,primerdf,correct){
    primerdf.new <- primerdf
    df.new <- df
    primerdf.new %<>% group_by(primerSet,species) %>% 
        mutate(maxEfficiency=max(efficiency,na.rm=TRUE)) %>% 
        ungroup() %>% 
        distinct(primerSet,species,maxEfficiency) %>% 
        group_by(primerSet) %>%
        mutate(maxByPrimer=max(maxEfficiency,na.rm=TRUE),medianByPrimer=median(maxEfficiency,na.rm=TRUE)) %>% 
        ungroup() %>% 
        mutate(corFactor=maxByPrimer/maxEfficiency,medianCorFactor=maxByPrimer/medianByPrimer) %>%
        mutate(primerSpp=paste(primerSet,species))
    df.new %<>% mutate(medianByPrimer=pull(primerdf.new,medianByPrimer)[match(primerSet,pull(primerdf.new,primerSet))],medianCorFactor=pull(primerdf.new,medianCorFactor)[match(primerSet,pull(primerdf.new,primerSet))]) %>%
        mutate(primerSpp=paste(primerSet,species)) %>% 
        mutate(maxEfficiency=pull(primerdf.new,maxEfficiency)[match(primerSpp,pull(primerdf.new,primerSpp))]) %>%
        mutate(maxEfficiency=if_else(is.na(maxEfficiency),medianByPrimer,maxEfficiency)) %>%
        mutate(corFactor=pull(primerdf.new,corFactor)[match(primerSpp,pull(primerdf.new,primerSpp))]) %>%
        mutate(corFactor=if_else(is.na(corFactor),medianCorFactor,corFactor)) #%>%
        #replace_na(list(corFactor=1)) 
        #select(-primerSpp,-meanEfficiency)
    #print(primerdf.new,n=Inf)
    if(correct==FALSE) {
        df.new <- df.new %>% select(-medianByPrimer,-medianCorFactor,-primerSpp,-corFactor)
    } else if(correct==TRUE) {
        df.new <- df.new %>% mutate(nreads=ceiling(nreads*corFactor)) %>% select(-medianByPrimer,-medianCorFactor,-primerSpp,-corFactor)
    } else stop(writeLines("correct must be TRUE or FALSE."))
    return(df.new)
}


# function to add a p value and r squared to a plot 
extract_p <- function(df,y,x,type,dp) {
    if(type=="lm") {
    fit <- df %>% lm(formula=paste(y,x,sep=" ~ "),data=.,na.action=na.omit) %>%
        broom::glance() %>% 
        mutate(pretty=paste0("R^2=",round(r.squared,digits=2),", ","p=",round(p.value,digits=dp))) %>%
        pull(pretty)
    } else if(type=="glm") {
    fit <- df %>% glm(formula=paste(y,x,sep=" ~ "),data=.,family=binomial(link="logit"),na.action=na.omit) %>%
        broom::tidy() %>% filter(term=="nReads") %>%
        mutate(pretty=paste0("p=",round(p.value,digits=dp))) %>%
        pull(pretty)
    } else stop(writeLines("type must be 'lm' or 'glm'."))
return(fit)
}


# make an NMDS plotting FUN
plot_nmds <- function(df,var) {
    n.var <- df %>% distinct(!!as.name(var)) %>% pull(!!as.name(var)) %>% length()
    if(n.var > 12) {
        cols <- colorRampPalette(ptol_pal()(12))(n.var)
    } else {
        cols <- ptol_pal()(n.var)
    }
    p <- df %>% ggplot(aes(x=NMDS1,y=NMDS2,color=!!as.name(var),fill=!!as.name(var))) + 
        geom_point(size=5,alpha=1,shape=24) +
        geom_convexhull(aes(fill=!!as.name(var),color=!!as.name(var)),alpha=0.1) +
        scale_discrete_manual(values=cols,aesthetic=c("color","fill")) + 
        theme_bw()
return(p)
}
