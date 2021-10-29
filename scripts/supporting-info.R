#!/usr/bin/env Rscript

# load libs and funs
source(here::here("scripts/load-data.R"))
library("xtable")


### GENERAL STATS ADULTS ###

# collapse eDNA
edna.collapsed <- edna.filt %>% 
    collapse_taxonomy(rmfw=TRUE,lifestage="Adult",collapse=FALSE) %>% 
    filter(partnerID=="MBA" & !grepl("WHIT",eventID)) %>% 
    filter(nreads>0)

# total reads & samples
edna.collapsed %>% pull(nreads) %>% sum()
edna.collapsed %>% distinct(sampleHash)

# get time duration
dur.edna <- edna.collapsed %>% distinct(eventDate) %>% arrange(eventDate) %>% pull(eventDate)

# filter trad 
trad.filtered <- trad.master %>% 
    collapse_taxonomy(rmfw=TRUE,lifestage="Adult",collapse=FALSE) %>%
    filter(partnerID=="MBA" & !grepl("WHIT",eventID))

# collapse taxonomy
trad.collapsed <- trad.filtered %>% 
    filter(grepl(" ",species)) %>%
    filter(eventDate %within% lubridate::interval(start=first(dur.edna),end=last(dur.edna)))

# get spp lists
trad.spp <- trad.collapsed %>% distinct(species) %>% pull()
edna.spp <- edna.collapsed %>% distinct(species) %>% pull()

# get trad numbers before time filter (seasonal PCoAs)
trad.filtered %>% distinct(species) %>% pull() %>% length()
trad.filtered %>% pull(individualCount) %>% sum()
trad.filtered %>% distinct(eventID,eventDate,fieldNumber) %>% count() %>% pull()

# get numbers after time duration filter (correlations models)
length(trad.spp)
trad.collapsed %>% pull(individualCount) %>% sum()
trad.collapsed %>% distinct(eventID,eventDate,fieldNumber) %>% count() %>% pull()
#
length(edna.spp)
edna.collapsed %>% pull(nreads) %>% sum()

# in trad, not eDNA
setdiff(trad.spp,edna.spp) %>% length()
# in eDNA not trad
setdiff(edna.spp,trad.spp) %>% length()

# shared spp
length(intersect(trad.spp,edna.spp))

# shared proportion
length(intersect(trad.spp,edna.spp)) / length(trad.spp)

# combined spp
length(union(trad.spp,edna.spp))

# species list slopes plot
trad.list <- trad.collapsed %>% group_by(species) %>% summarise(totalTrad=sum(individualCount),.groups="drop")
edna.list <- edna.collapsed %>% group_by(species) %>% summarise(totalEdna=sum(nreads),.groups="drop")

# join and clean
lists.joined <- full_join(trad.list,edna.list) %>% 
    mutate_all(~replace_na(.,0))


# make the ranking
spp.ranked <- lists.joined %>% 
    mutate(trad=row_number(desc(totalTrad)),edna=row_number(desc(totalEdna))) %>% 
    select(species,edna,trad) %>% 
    arrange(edna) %>% 
    mutate(rankDiff=edna-trad) %>%
    pivot_longer(cols=c(-species,-rankDiff),names_to="data",values_to="rank")

# plot
p <- spp.ranked %>% ggplot(aes(x=data,y=rank,group=species)) + 
    geom_line(aes(color=rankDiff),size=1.3,alpha=0.75) + 
    geom_text(data=filter(spp.ranked,data=="edna"),aes(label=species),check_overlap=TRUE,hjust=1,nudge_x=-0.02) +
    geom_text(data=filter(spp.ranked,data=="trad"),aes(label=species),check_overlap=TRUE,hjust=0,nudge_x=0.02) +
    scale_y_reverse() +
    theme_void() +
    theme(legend.position="none") +
    scale_colour_distiller(palette="RdYlBu")
plot(p)
ggsave(filename=here("temp/results/figures/species-slopes.svg"),plot=p,width=250,height=500,units="mm")

# check the zero spp
lists.joined %>% filter(totalEdna==0) %>% arrange(species)
lists.joined %>% filter(totalTrad==0) %>% arrange(species) %>% print(n=Inf)

# annotate and print out the list
lists.joined %>% 
    left_join(rename(tax.table,species=sciNameValid)) %>% 
    mutate(speciesUncollapsed=species) %>% 
    adult_edit(collapse=TRUE) %>% 
    mutate(species=if_else(species==speciesUncollapsed,"",species)) %>%
    rename(collapse=species,species=speciesUncollapsed) %>%
    arrange(class,order,family,species) %>% 
    select(family,species,collapse,totalEdna,totalTrad) %>%
    mutate(species=str_replace_all(species,"^","\\\\emph{"),species=str_replace_all(species,"$","}")) %>%
    rename(Family=family,Species=species,`Merged taxa as`=collapse,`Total eDNA reads`=totalEdna,`Total trawl abundance`=totalTrad) %>%
    mutate_if(is.numeric,~prettyNum(.,big.mark=",")) %>%
    xtable(caption="blahhh", digits=c(0,0,0,0,0,0)) %>%
    print.xtable(include.rownames=FALSE,booktabs=TRUE,sanitize.text.function=identity,caption.placement="top",size="tiny")

# sum up
lists.joined %>% summarise(totalTrad=sum(totalTrad),totalEdna=sum(totalEdna))


### BIOINFORMATICS SUMMARY TABLE ###

# load and format
read_csv(here("meta-fish-pipe/results/stats-summary.csv")) %>% 
    pivot_wider(names_from=library,values_from=nreads) %>%
    mutate_all(~prettyNum(.,big.mark=",")) %>%
    xtable(caption="blahhh", digits=c(0,0,0,0,0,0)) %>%
    print.xtable(include.rownames=FALSE,booktabs=TRUE,sanitize.text.function=identity,caption.placement="top")


### SEQUENCING DEPTH ###

# plot 
p <- edna.filt %>% filter(primerSet=="tele02" & partnerID=="MBA") %>%
    mutate(month=month(eventDate,label=TRUE),year=year(eventDate),yearMonth=paste(year,month,sep="-")) %>%
    mutate(yearMonth=fct_relevel(yearMonth,c("2017-Feb","2017-Mar","2017-Apr","2017-May","2017-Jun","2017-Jul","2017-Aug","2017-Sep","2017-Oct","2017-Nov","2018-Feb","2018-Mar","2018-Apr"))) %>%
    group_by(library,sampleHash,yearMonth) %>% 
    summarise(nreads=sum(nreads),.groups="drop") %>%
        ggplot(aes(x=yearMonth,y=nreads)) +
        geom_boxplot(outlier.shape=NA,color="darkgrey") +
        geom_jitter(aes(color=yearMonth),width=0.3,height=0) +# ,position="jitter"
        scale_y_log10() +
        scale_x_discrete(labels=c("2017\nFeb","2017\nMar","2017\nApr","2017\nMay","2017\nJun","2017\nJul","2017\nAug","2017\nSep","2017\nOct","2017\nNov","2018\nFeb","2018\nMar","2018\nApr")) +
        ggthemes::theme_clean(base_size=12) +
        scale_color_manual(values=colorRampPalette(ptol_pal()(12))(13)) +
        labs(x="Month",y="Read depth (log10)") +
        theme(legend.position="none")
plot(p)
ggsave(filename=here("temp/results/figures/read-depth-libs.svg"),plot=p,width=200,height=105,units="mm")


### NMDS BY DEPTH/LOCATION ###
#library("vegan")
#renv::install("cmartin/ggConvexHull")
library("ggConvexHull")

# collapse
trad.collapsed <- trad.master %>% collapse_taxonomy(rmfw=TRUE,lifestage="Adult",collapse=FALSE)# ORDINATIONS
edna.collapsed <- edna.filt %>% collapse_taxonomy(rmfw=TRUE,lifestage="Adult",collapse=FALSE)# ORDINATIONS

# make a list of species common to edna and trad (by partnerID)
comb.spp <- full_join(distinct(edna.collapsed,partnerID,species),distinct(trad.collapsed,partnerID,species)) %>% arrange(partnerID,species)

# expand for all spp and get groups
edna.expanded <- edna.collapsed %>% expand_and_summarise(sppdf=comb.spp,eventdf=events,primerdf=primer.bias,method="edna",spatialgroup="locality",temporalgroup="day",correct=FALSE)# EGGS/ADULTS/SPAWNING

# format df
edna.org <- edna.expanded %>% 
    filter(partnerID=="MBA" & localityID!="WHIT") %>% 
    group_by(species) %>% 
    mutate(sppn=sum(nReads)) %>% 
    ungroup() %>%
    filter(sppn>0) %>%
    mutate(label=paste(localitySite,lubridate::month(eventDate,label=TRUE),eventDate,replicateFilter,library,sep=".")) 

# make matrix
edna.mat.mds <- edna.org %>%
    distinct(label,nReads,species) %>%
    pivot_wider(names_from=label,values_from=nReads) %>%
    data.frame(row.names=1) %>% 
    t() %>%
    as.matrix()

# run NMDS
edna.mds <- edna.mat.mds %>% metaMDS(autotransform=TRUE,k=2,distance="bray",binary=FALSE,trace=3,trymax=500)
edna.mds %>% print()

# extract, make factors
nmds.tbl <- edna.mds %>% 
    vegan::scores() %>%
    as_tibble(rownames="Sample") %>%
    mutate(yearMonth=paste(str_split_fixed(Sample,"\\.",8)[,4],str_split_fixed(Sample,"\\.",8)[,3],sep="-"),locality=str_split_fixed(Sample,"\\.",8)[,1],site=paste(str_split_fixed(Sample,"\\.",8)[,1],str_split_fixed(Sample,"\\.",8)[,2],sep="-"),depth=str_split_fixed(Sample,"\\.",8)[,2],lib=str_split_fixed(Sample,"\\.",8)[,8]) %>%
    mutate(yearMonth=fct_relevel(yearMonth,c("2017-Feb","2017-Mar","2017-Apr","2017-May","2017-Jun","2017-Jul","2017-Aug","2017-Sep","2017-Oct","2017-Nov","2018-Feb","2018-Mar","2018-Apr"))) 

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

# plot all by YEARMONTH 
p <- nmds.tbl %>% plot_nmds(var="yearMonth")
ggsave(filename=here("temp/results/figures/months-nmds.svg"),plot=p,width=200,height=170,units="mm")


# plot all by LIBRARY
p <- nmds.tbl %>% plot_nmds(var="lib")
ggsave(filename=here("temp/results/figures/libs-nmds.svg"),plot=p,width=200,height=200,units="mm")


# plot all by LOCALITY
p <- nmds.tbl %>% plot_nmds(var="locality")
ggsave(filename=here("temp/results/figures/locality-nmds.svg"),plot=p,width=200,height=170,units="mm")


# plot all by SITE
p <- nmds.tbl %>% plot_nmds(var="site")
ggsave(filename=here("temp/results/figures/site-nmds.svg"),plot=p,width=200,height=170,units="mm")


# plot all by DEPTH
p <- nmds.tbl %>% plot_nmds(var="depth")
ggsave(filename=here("temp/results/figures/depth-nmds.svg"),plot=p,width=200,height=170,units="mm")


### NEGATIVE CONTROLS ###

# table of negative controls
read_csv(here("meta-fish-pipe/results/controls-summary.csv")) %>% 
    mutate(partnerID=str_split_fixed(eventID,"-",4)[,1]) %>% 
    filter(partnerID!="EA") %>%
    group_by(library,blankType,sampleHash) %>%
    summarise(sampleReads=sum(nReads),.groups="drop") %>%
    arrange(library,blankType,desc(sampleReads)) %>% 
    mutate(blankType=str_replace_all(blankType,"Blank"," blank"), blankType=str_replace_all(blankType,"Well","Tag")) %>%
    rename(Library=library,`Control type`=blankType,`Sample ID`=sampleHash,`Number reads`=sampleReads) %>%
    mutate_if(is.numeric,~prettyNum(.,big.mark=",")) %>%
    xtable(caption="blahhh") %>% 
    print.xtable(include.rownames=FALSE,booktabs=TRUE,sanitize.text.function=identity,caption.placement="top",size="scriptsize")

### READS AND SAMPLES ###

# table of samples
edna.filt %>% filter(partnerID=="MBA") %>%
    group_by(library,eventID,sampleHash) %>%
    summarise(sampleReads=sum(nreads),.groups="drop") %>%
    group_by(library) %>% 
    summarise(nSamples=n_distinct(sampleHash),meanReads=mean(sampleReads),sd=sd(sampleReads),.groups="drop") %>%
    mutate(meanReads=round(meanReads,digits=0),sd=round(sd,digits=0)) %>%
    rename(Library=library,`Samples (n)`=nSamples,`Reads (mean)`=meanReads,`Reads (sd)`=sd) %>%
    mutate_if(is.numeric,~prettyNum(.,big.mark=",")) %>%
    xtable(caption="blahhh") %>% 
    print.xtable(include.rownames=FALSE,booktabs=TRUE,sanitize.text.function=identity,caption.placement="top")
