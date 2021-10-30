#!/usr/bin/env Rscript

# load libs
suppressPackageStartupMessages({
    library("here")
    library("tidyverse")
    library("lubridate")
    #library("RColorBrewer")
    library("magrittr")
    library("ggthemes")
    library("FactoMineR")
    #library("ggridges")
    library("broom")
    library("rfishbase")
    library("glmmTMB")
    library("DHARMa")
    library("emmeans")
    library("performance")
    library("vegan")
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

# 
# fun to edit the species for traits analyses
#traits_edit <- function(df) {
#    df.ed <- df
#return(df.ed)
#}


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


## function to process and collapse the trad and eDNA data
#process_annotate <- function(df,eventdf,method,summary_func,spatialgroup,temporalgroup){
#    events.sub <- eventdf %>% distinct(eventID,localityID,localitySite,year,month,day)
#    if(method=="traditional") {
#        completed.df <- df %>% 
#            group_by(partnerID,eventID,eventDate,fieldNumber,species) %>% 
#            summarise(individualCount=sum(individualCount),weightInGrams=sum(weightInGrams),.groups="drop") %>% # merge dup species after rename
#            group_by(partnerID) %>%
#            tidyr::complete(tidyr::nesting(eventID,eventDate,fieldNumber),species,fill=list(individualCount=0,weightInGrams=0)) %>% # add missing zeros for all species by partnerID. Nesting allows only those combinations of samples present in the dataset
#            ungroup() %>%
#            mutate(localityID=pull(events.sub,localityID)[match(eventID,pull(events.sub,eventID))],# add the localitySite details
#            localitySite=pull(events.sub,localitySite)[match(eventID,pull(events.sub,eventID))],
#            year=pull(events.sub,year)[match(eventID,pull(events.sub,eventID))],
#            month=pull(events.sub,month)[match(eventID,pull(events.sub,eventID))],
#            day=pull(events.sub,day)[match(eventID,pull(events.sub,eventID))]) 
#        #
#        if(temporalgroup=="month") {
#            completed.df %<>% mutate(temporalGroup=paste(year,month,"1",sep="-")) } 
#        else if(temporalgroup=="day"){
#            completed.df %<>% mutate(temporalGroup=paste(year,month,day,sep="-")) } 
#        else stop(writeLines("temporalgroup must be 'month' or 'day'."))
#        if(spatialgroup=="locality") {
#            completed.df %<>% mutate(spatialGroup=localityID) } 
#        else if(spatialgroup=="partner") {
#            completed.df %<>% mutate(spatialGroup=partnerID) } 
#        else if(spatialgroup=="site"){
#            completed.df %<>% mutate(spatialGroup=localitySite) } 
#        else stop(writeLines("spatialgroup must be 'locality', 'site', or 'partner'."))
#        #
#        completed.df %<>%
#            group_by(partnerID,spatialGroup,temporalGroup,species) %>% # collapse to get summarised numbers/grams per location period
#            summarise(individualsByGroup=ceiling(summary_func(individualCount,na.rm=TRUE)),weightInGramsByGroup=ceiling(summary_func(weightInGrams,na.rm=TRUE)),nReplicates=n(),.groups="drop") %>% 
#            group_by(partnerID,spatialGroup,temporalGroup) %>%  #
#            mutate(individualsByGroupTotal=sum(individualsByGroup,na.rm=TRUE),weightInGramsByGroupTotal=sum(weightInGramsByGroup)) %>% # get totals over events (after mean)
#            ungroup() %>% 
#            mutate(individualsByGroupProportion=individualsByGroup/individualsByGroupTotal,individualsByGroupRate=individualsByGroup/nReplicates) %>% # get prop INDIVS, prop WEIGHTS
#            mutate(weightInGramsByGroupProportion=weightInGramsByGroup/weightInGramsByGroupTotal,weightInGramsByGroupRate=weightInGramsByGroup/nReplicates) %>%
#            mutate(weightInGramsByGroup=if_else(partnerID=="EA",NaN,weightInGramsByGroup),weightInGramsByGroupTotal=if_else(partnerID=="EA",NaN,weightInGramsByGroupTotal),weightInGramsByGroupRate=if_else(partnerID=="EA",NaN,weightInGramsByGroupRate)) %>% # clean nans
#            mutate_all(~replace(.,is.nan(.),NA)) %>% # clean nans
#            #select(partnerID,location,period,species,individualCount,eventIndividualCount,eventIndividualProportion,weightInGrams,eventWeightInGrams,eventWeightInGramsProportion,nReplicates) %>% 
#            mutate(methodology="traditional")
#    #
#    } else if(method=="edna") {
#        completed.df <- df %>% 
#            group_by(primerSet,partnerID,eventID,eventDate,replicateFilter,replicatePCR,sampleHash,species) %>% 
#            summarise(nreads=sum(nreads,na.rm=TRUE),.groups="drop") %>% # merge dup species after rename
#            group_by(partnerID) %>%
#            tidyr::complete(tidyr::nesting(primerSet,eventID,eventDate,replicateFilter,replicatePCR,sampleHash),species,fill=list(nreads=0)) %>% # add missing zeros for all species by partnerID. Nesting allows only those combinations of samples present in the dataset
#            ungroup() %>% 
#            mutate(localityID=pull(events.sub,localityID)[match(eventID,pull(events.sub,eventID))],# add the localitySite details
#            localitySite=pull(events.sub,localitySite)[match(eventID,pull(events.sub,eventID))],
#            year=pull(events.sub,year)[match(eventID,pull(events.sub,eventID))],
#            month=pull(events.sub,month)[match(eventID,pull(events.sub,eventID))],
#            day=pull(events.sub,day)[match(eventID,pull(events.sub,eventID))])
#        #
#        if(temporalgroup=="month") {
#            completed.df %<>% mutate(temporalGroup=paste(year,month,"1",sep="-")) } 
#        else if(temporalgroup=="day"){
#            completed.df %<>% mutate(temporalGroup=paste(year,month,day,sep="-")) } 
#        else stop(writeLines("temporalgroup must be 'month' or 'day'."))
#        if(spatialgroup=="locality") {
#            completed.df %<>% mutate(spatialGroup=localityID) } 
#        else if(spatialgroup=="partner") {
#            completed.df %<>% mutate(spatialGroup=partnerID) } 
#        else if(spatialgroup=="site"){
#            completed.df %<>% mutate(spatialGroup=localitySite) } 
#        else stop(writeLines("spatialgroup must be 'locality', 'site', or 'partner'."))
#        #
#        completed.df %<>%
#            rename(nReads=nreads) %>%
#            group_by(sampleHash) %>%  # get sum of reads per sample replicate
#            mutate(readsBySampleTotal=sum(nReads,na.rm=TRUE)) %>% 
#            ungroup() %>% 
#            group_by(partnerID,spatialGroup,temporalGroup,species) %>% # collapse to get summarised collapsed reads per location/period
#            mutate(readsByGroup=ceiling(summary_func(nReads,na.rm=TRUE)),nGroupReplicates=n()) %>%
#            ungroup() %>% 
#            group_by(partnerID,spatialGroup,temporalGroup) %>%  # get sum of collapsed reads per event
#            mutate(readsByGroupTotal=sum(nReads,na.rm=TRUE)) %>% # get totals over events (after mean)
#            ungroup() %>% 
#            mutate(readsByGroupProportion=readsByGroup/readsByGroupTotal) %>% # get prop READS
#            mutate(methodology="edna") #%>%
#            #select(methodology,primerSet,partnerID,eventID,localityID,localitySite) %>% 
#    } else stop(writeLines("Method must be 'traditional' or 'edna'."))
#    return(completed.df)
#}


## function to join the eDNA and trad survey dataframes
#join_surveys <- function(trad,edna,comparison,minreads){
#    edna.red <- edna %>% rename(count=nReads,eventCount=eventReadCount,eventProportion=eventReadProportion) %>% filter(eventCount>=minreads)
#    if(comparison=="grams") { 
#        trad.red <- trad %>% rename(count=weightInGramsByGroup,eventCount=eventWeightInGrams,eventProportion=eventWeightInGramsProportion) %>% 
#            select(-individualCount,-eventIndividualCount,-eventIndividualProportion)
#        }
#    else if(comparison=="individuals") {
#        trad.red <- trad %>% rename(count=individualCount,eventCount=eventIndividualCount,eventProportion=eventIndividualProportion) %>%
#            select(-weightInGrams,-eventWeightInGrams,-eventWeightInGramsProportion)
#        }
#    else stop(writeLines("Comparison must be 'grams' or individuals'."))
#    trad.edna.comb <- bind_rows(edna.red,trad.red)
#    trad.edna.comb %<>% group_by(partnerID) %>% # now add all the species missing from either the trad or eDNA surveys ... makes some empty events 
#            tidyr::complete(tidyr::nesting(location,period,methodology),species,fill=list(count=0,eventProportion=0)) %>% # fill any remaining combinations not present
#            ungroup() %>% 
#            group_by(location,period,methodology) %>% 
#            fill(nReplicates,eventCount,.direction="downup") %>% ungroup() %>%
#            replace_na(list(eventCount=0,nReplicates=0))
#    return(trad.edna.comb)
#}


## function to rank species and subset the dataframe by date
#rank_subset <- function(df,time,rankmethod,n){
#    df.rank <- df %>% group_by(species) %>%
#        mutate(sppProps=sum(count[which(methodology==rankmethod)],na.rm=TRUE)) %>% ungroup() %>%
#        mutate(rank=dense_rank(desc(sppProps))) %>% # dense_rank(desc(sppProps)) # dense_rank(sppProps)
#        select(-sppProps) %>% 
#        mutate(period=ymd(period))
#    if(n=="all"){
#        df.rank <- df.rank
#        }
#    else {
#        df.rank %<>% filter(rank %in% 1:n) 
#        }
#    dur.edna <- df.rank %>% filter(methodology=="edna" & eventCount > 0) %>% distinct(period) %>% arrange(period) %>% pull(period)
#    if(time=="duration"){
#        df.sub <- df.rank %>% filter(period %within% lubridate::interval(start=first(dur.edna),end=last(dur.edna)))
#        }
#    else if(time=="identical"){
#        df.sub <- df.rank %>% filter(period %in% dur.edna)
#        }
#    else if(time=="all"){
#        df.sub <- df.rank
#        }
#    else stop(writeLines("Time must be: 'duration', 'identical', or 'all'."))
#return(df.sub)
#}


## function to plot a ridge plot
#plot_ridges <- function(df){
#    pl <- df %>% ggplot(aes(x=period,y=0, height=asin(sqrt(eventProportion)),fill=fct_relevel(methodology,"traditional","edna"))) +
#            geom_ridgeline(alpha=0.5,linetype="blank") +
#            facet_wrap(.~fct_reorder(species,rank),scales="free") +
#            theme_clean(base_size=10) +
#            scale_fill_ptol(name=NULL)
#return(pl)
#}


## function to plot correlations
#plot_corrs <- function(df,comparison,transform_func,partition){
#    if(partition==TRUE) {
#        p <- df %>% ggplot(aes(y=transform_func(!!as.name(paste0(comparison,".edna"))),x=transform_func(!!as.name(paste0(comparison,".traditional"))),color=location)) +
#            geom_smooth(method="lm",se=FALSE) +
#            geom_point() +
#            facet_wrap(.~fct_reorder(species,rank),scales="free") +
#            labs(x="Proportion of traditional survey catch (arcsine squareroot transformed)",y="Proportion of eDNA survey catch (arcsine squareroot transformed)") +
#            theme_clean(base_size=10) +
#            scale_color_ptol()
#    } else if (partition==FALSE) {
#        p <- df %>% ggplot(aes(y=transform_func(!!as.name(paste0(comparison,".edna"))),x=transform_func(!!as.name(paste0(comparison,".traditional"))))) +
#            geom_smooth(method="lm",se=TRUE,color="#4477AA",fill="grey90") +
#            geom_point() +
#            facet_wrap(.~fct_reorder(species,rank),scales="free") +
#            labs(x="Proportion of traditional survey catch (arcsine squareroot transformed)",y="Proportion of eDNA survey catch (arcsine squareroot transformed)") +
#            theme_clean(base_size=10)
#    } else stop(writeLines("partition must be: 'TRUE' or 'FALSE'."))
#    return(p)
#}


## function to plot glms
#plot_glms <- function(df,compy,compx,transform_func_x,partition,label_y,label_x){
#    df <- df %>% drop_na(any_of(compy)) %>% mutate(binvar=if_else(!!as.name(compy)==TRUE,1,0))
#    if(partition==TRUE) {
#        p <- df %>% 
#        ggplot(aes(x=transform_func_x(!!as.name(compx)),y=binvar)) +
#            geom_smooth(method="glm",method.args=list(family=binomial(link="logit"),na.action=na.omit),se=TRUE,alpha=0.5,color="#4477AA") +
#            geom_rug(data=filter(df,binvar==0),sides="b",alpha=1) +
#            geom_rug(data=filter(df,binvar==1),sides="t",alpha=1) +
#            #geom_line(aes(y=.fitted),color="#4477AA",size=0.75) +
#            #geom_ribbon(aes(ymin=.fitted-1.96*.se.fit, ymax=.fitted+1.96*.se.fit),alpha=0.1) +
#            facet_wrap(.~fct_reorder(species,rank),scales="free") +
#            labs(x=label_x,y=label_y) + # x="Proportion of total eDNA event reads (arcsine squareroot transformed)",y="Pr (breeding month)"
#            ggthemes::theme_clean(base_size=8)
#    } else if (partition==FALSE) {
#        p <- df %>% 
#        ggplot(aes(x=transform_func_x(!!as.name(compx)),y=binvar)) +
#            geom_smooth(method="glm",method.args=list(family=binomial(link="logit"),na.action=na.omit),se=TRUE,alpha=0.2,color="#4477AA") +
#            geom_rug(data=filter(df,binvar==0),sides="b",alpha=0.2) +
#            geom_rug(data=filter(df,binvar==1),sides="t",alpha=0.2) +
#            #geom_line(aes(y=.fitted),color="#4477AA",size=0.75) +
#            #geom_ribbon(aes(ymin=.fitted-1.96*.se.fit, ymax=.fitted+1.96*.se.fit),alpha=0.1) +
#            labs(x=label_x,y=label_y) +
#            ggthemes::theme_clean(base_size=16)
#    } else stop(writeLines("partition must be: 'TRUE' or 'FALSE'."))
#    return(p)
#}


# function to plot boxes for breeding
#plot_trait_boxes <- function(df,compy,compx,transform_func_y,partition,label_y,label_x){
#    df <- df %>% drop_na(any_of(compx))
#    if(partition==TRUE) {
#        p <- df %>% ggplot(aes(x=!!as.name(compx),y=transform_func_y(!!as.name(compy)))) +
#            #geom_smooth(method="lm",se=TRUE,color="#4477AA") +
#            geom_boxplot() +
#            facet_wrap(.~fct_reorder(species,rank),scales="free") +
#            labs(x=label_x,y=label_y) +
#            ggthemes::theme_clean(base_size=8)
#    } else if (partition==FALSE) {
#        p <- df %>% ggplot(aes(x=!!as.name(compx),y=transform_func_y(!!as.name(compy)))) +
#            #geom_smooth(method="lm",se=TRUE,color="#4477AA") +
#            geom_boxplot() +
#            labs(x=label_x,y=label_y) +
#            ggthemes::theme_clean(base_size=16)
#    } else stop(writeLines("partition must be: 'TRUE' or 'FALSE'."))
#    return(p)
#}


# function to plot breeding correlations
#plot_trait_lms <- function(df,compy,compx,transform_func_y,transform_func_x,label_y,label_x,partition){
#    df <- df %>% drop_na(any_of(compx))
#    if(partition==TRUE) {
#        p <- df %>% ggplot(aes(x=transform_func_x(!!as.name(compx)),y=transform_func_y(!!as.name(compy)))) +
#            geom_smooth(method="lm",se=TRUE,color="#4477AA",na.rm=TRUE) +
#            geom_point(na.rm=TRUE) +
#            facet_wrap(.~fct_reorder(species,rank),scales="free") +
#            labs(x=label_x,y=label_y) +
#            ggthemes::theme_clean(base_size=8)
#    } else if (partition==FALSE) {
#        p <- df %>% ggplot(aes(x=transform_func_x(!!as.name(compx)),y=transform_func_y(!!as.name(compy)))) +
#            geom_smooth(method="lm",se=TRUE,color="#4477AA",na.rm=TRUE) +
#            geom_point(na.rm=TRUE) +
#            labs(x=label_x,y=label_y) +
#            ggthemes::theme_clean(base_size=16)
#    } else stop(writeLines("partition must be: 'TRUE' or 'FALSE'."))
#    return(p)
#}


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


## arcsine squareroot fun
#asin_sqrt <- function(x){
#    xt <- asin(sqrt(x))
#    return(xt)
#}


# log plus1 fun
logplus1 <- function(x){
    xt <- log(x+1)
    return(xt)
}


## do nothing fun
#do_nothing <- function(x){
#    xt <- x
#    return(xt)
#}


# 4th root function
fourth_root <- function(x){
    xt <- x^(1/4)
    return(xt)
}


## function to extract regression params
#regression_extract <- function(df,summary_func,transform_func,comparison){
#    df.spp <- df %>% 
#        group_by(species) %>% 
#        mutate(summary.traditional=summary_func(count.traditional,na.rm=TRUE),summary.edna=summary_func(count.edna,na.rm=TRUE)) %>%
#        tidyr::nest() %>% 
#        mutate(fit=map(data, ~ lm(transform_func(!!as.name(paste0(comparison,".edna"))) ~ transform_func(!!as.name(paste0(comparison,".traditional"))), na.action="na.exclude", data = .x)), 
#            tidied = map(fit,broom::tidy),
#            glanced = map(fit,broom::glance)) %>%
#        ungroup() %>% 
#        tidyr::unnest(glanced) %>%
#        select(-p.value, -statistic) %>%
#        tidyr::unnest(tidied) %>%
#        filter(term!="(Intercept)") %>% 
#        tidyr::unnest(data) %>%
#        distinct(species,rank,estimate,p.value,r.squared,summary.traditional,summary.edna) %>%
#        arrange(rank) %>%
#        mutate(significant=if_else(p.value<=0.05,TRUE,FALSE)) %>%
#        mutate_if(is.numeric, ~replace(., is.nan(.), NA))
#        # 
#    df.t <- df %>% mutate(comp.edna.trans=transform_func(!!as.name(paste0(comparison,".edna"))), comp.trad.trans=transform_func(!!as.name(paste0(comparison,".traditional"))))
#    df.aug <- lm(formula=comp.edna.trans~comp.trad.trans,na.action="na.exclude",data=df.t) %>% 
#        augment(data=df.t) %>% 
#        group_by(species) %>%  
#        summarise(mean.resid=mean(.resid,na.rm=TRUE),se.resid=se(.resid))
#    df.comb <- left_join(df.spp,df.aug)
#    return(df.comb)
#}


# standard error function
se <- function(x) sqrt(var(x)/length(x))


## function to add fishbase data to table
#add_fishbase <- function(df) {
#    df.sp <- df %>% select(species) %>% 
#        mutate(speciesFb=species) %>%
#        mutate(speciesFb=str_replace_all(speciesFb,"Merlangius/Melanogrammus","Merlangius merlangus"), 
#            speciesFb=str_replace_all(speciesFb,"Hyperoplus/Ammodytes","Ammodytes tobianus"),
#            speciesFb=str_replace_all(speciesFb,"Chelon labrosus/ramada","Chelon ramada"),
#            speciesFb=str_replace_all(speciesFb,"Triglidae","Eutrigla gurnardus"), 
#            speciesFb=str_replace_all(speciesFb,"Pomatoschistus microps/minutus","Pomatoschistus minutus"),
#            speciesFb=str_replace_all(speciesFb,"Aphia/Crystallogobius","Aphia minuta"))
#    fb.df <- rfishbase::species(pull(df.sp,speciesFb),server="fishbase") %>% 
#        select(Species,DemersPelag) %>% 
#        mutate(DemersPelag=str_replace_all(DemersPelag,"bathydemersal","demersal")) %>%
#        rename(speciesFb=Species) %>% 
#        left_join(df.sp) %>% 
#        select(-speciesFb) %>% 
#        rename(lifestyle=DemersPelag)
#    df.ann <- left_join(df,fb.df)
#    return(df.ann)
#}


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
    custom.tab <- read_csv(file=here("assets/fishbase-spawning-custom.csv")) %>% 
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
        mutate(fecundityMinMean=mean(fecundityMin,na.rm=TRUE),fecundityMaxMean=mean(fecundityMax,na.rm=TRUE),fecundityMean=mean(c(fecundityMinMean,fecundityMaxMean),na.rm=TRUE)) %>%
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
    df.joined <- df %>% left_join(spawn.tab) %>% left_join(custom.tab) %>% left_join(larvae.tab) %>% left_join(fec.tab) %>% left_join(life.tab) %>% left_join(length.tab) %>% left_join(ecology.tab)
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
    trad.edna.comb <- left_join(edna,trad) %>% 
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
            full_join(sppdf) %>% # add the names from trad surveys
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


## a pseudo r2 function
#r2.corr.mer <- function(m) {
#   lmfit <-  lm(model.response(model.frame(m)) ~ fitted(m))
#   summary(lmfit)$r.squared
#}


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
