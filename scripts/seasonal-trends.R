#!/usr/bin/env Rscript

# load libs and funs
source(here::here("scripts/load-data.R"))

# make a figures output directory if needed
if(!dir.exists(here("temp/results/figures"))){
    dir.create(here("temp/results/figures"),recursive=TRUE)
}

#### EDNA ####

# collapse
edna.collapsed <- edna.filt %>% collapse_taxonomy(rmfw=TRUE,lifestage="Adult",collapse=FALSE)

# expand for all spp and get groups
edna.expanded <- edna.collapsed %>% expand_and_summarise(sppdf=distinct(edna.collapsed,species),eventdf=events,primerdf=primer.bias,method="edna",spatialgroup="locality",temporalgroup="day",correct=TRUE)# EGGS/ADULTS/SPAWNING

# filter format df
edna.org <- edna.expanded %>% 
    group_by(species) %>% 
    mutate(sppn=sum(nReads)) %>% 
    ungroup() %>%
    filter(sppn>1000) %>% 
    mutate(label=paste(localitySite,lubridate::month(eventDate,label=TRUE),eventDate,replicateFilter,sep=".")) %>%
    distinct(label,species,nReads,readsBySampleTotal)

# get stats
edna.org %>% distinct(species) %>% print(n=Inf)
edna.org %>% distinct(label) %>% print(n=Inf)
edna.org %>% mutate(date=lubridate::ymd(str_split_fixed(label,"\\.",4)[,3])) %>% distinct(date) %>% arrange(date) %>% print(n=Inf)


# run linear models for all
edna.org.lm <- edna.org %>% 
        group_by(species) %>% 
        tidyr::nest() %>% 
        mutate(fit=purrr::map(data, ~ lm(logplus1(nReads)~logplus1(readsBySampleTotal),na.action="na.exclude",data=.x)), 
            aug=purrr::map(fit,broom::augment)) %>%
        ungroup() %>% 
        tidyr::unnest(c(data,aug)) %>% 
        mutate(nReads=logplus1(nReads)) %>% 
        distinct(label,species,nReads,.std.resid) %>%
        arrange(label,species)

# make matrix (log reads)
edna.mat <- edna.org.lm %>% 
    distinct(label,species,nReads) %>%
    pivot_wider(names_from=label,values_from=nReads) %>%
    data.frame(row.names=1) %>% 
    t() %>%
    as.matrix()

# PCoA - Principal coordinates analysis
edna.mat.dist <- vegdist(edna.mat,method="euclidean")
edna.mat.pcoa <- cmdscale(edna.mat.dist)

# convert and format (pcoa)
edna.pc1.tib <- tibble(label=rownames(edna.mat.pcoa), PC1=edna.mat.pcoa[,2]) %>% # USE PC2
    mutate(yearMonth=str_split_fixed(label,"\\.",6)[,3]) %>%
    mutate(yearMonth=fct_relevel(yearMonth,c("Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov"))) 


# plot
p <- edna.pc1.tib %>%
        ggplot(aes(x=yearMonth,y=PC1)) +
        geom_point(color="gray30",size=2.5,alpha=0.8,shape=16) + 
        geom_smooth(aes(x=as.numeric(yearMonth)),method="gam",formula=y~s(x,bs="cs",k=4),color="#737145",alpha=0.5,fill="gray90") + 
        ggthemes::theme_clean(base_size=12) +
        labs(x="Date", y="eDNA PCoA-2 scores")
#plot(p)
ggsave(filename=here("temp/results/figures/edna-pcoa.svg"),plot=p,width=130,height=80,units="mm")

# SUPPORTING - plot edna PCA VARIABLES 
p <- plot.PCA(PCA(edna.mat,graph=FALSE,scale.unit=FALSE),choix="var",select="contrib 12",ggoptions=list(size=2,circle.lwd=0.3,circle.lty=3,line.lty=3,line.lwd=0.3)) + theme(plot.title=element_blank(),axis.title=element_text(size=8))
#plot(p)
ggsave(filename=here("temp/results/figures/edna-pcoa-plot-resids.svg"),plot=p,width=130,height=80,units="mm")


##### TRAWL #####

# filter mba, non-species, fw, and uneeded cols
trad.collapsed <- trad.master %>% filter(partnerID=="MBA" & grepl("L4",eventID)) %>%
    collapse_taxonomy(rmfw=TRUE,lifestage="Adult",collapse=FALSE) %>% 
    select(-basisOfRecord,-partnerID,-lifestage)

# create haul labels fill, and get tots per haul
trad.org <- trad.collapsed %>% 
    group_by(eventID) %>%
    mutate(haul=dense_rank(fieldNumber)) %>% 
    ungroup() %>% 
    mutate(haul=str_replace_all(haul,"1","A"),haul=str_replace_all(haul,"2","B"),haul=str_replace_all(haul,"3","C")) %>% 
    mutate(label=paste("L4-trawl",lubridate::month(eventDate,label=TRUE),eventDate,haul,sep=".")) %>% 
    select(label,species,individualCount,weightInGrams) %>%
    tidyr::complete(tidyr::nesting(label),species,fill=list(individualCount=0,weightInGrams=0)) %>% 
    group_by(label) %>%
    mutate(individualCountByHaul=sum(individualCount),weightInGramsByHaul=sum(weightInGrams)) %>%
    ungroup() %>%
    arrange(label,species)


# get stats
trad.org %>% distinct(species) %>% print(n=Inf)
trad.org %>% distinct(label) %>% print(n=Inf)
trad.org %>% mutate(date=lubridate::ymd(str_split_fixed(label,"\\.",4)[,3])) %>% distinct(date) %>% arrange(date) %>% print(n=Inf)
trad.org %>% summarise(sum=sum(individualCount))

# make a lm per spp with log+1 counts
trad.org.lm <- trad.org %>% 
        group_by(species) %>% 
        tidyr::nest() %>% 
        mutate(fit=purrr::map(data, ~ lm(logplus1(individualCount)~logplus1(individualCountByHaul),na.action="na.exclude",data=.x)), 
            aug=purrr::map(fit,broom::augment)) %>%
        ungroup() %>% 
        tidyr::unnest(c(data,aug)) %>% 
        select(species,label,individualCount,weightInGrams,.std.resid) %>%
        mutate(individualCount=logplus1(individualCount),weightInGrams=logplus1(weightInGrams))

# make matrix for LOG INDIVS
trad.mat <- trad.org.lm %>% 
    distinct(label,individualCount,species) %>%
    pivot_wider(names_from=label,values_from=individualCount) %>%
    data.frame(row.names=1) %>% 
    t() %>%
    as.matrix()

# try PCoA
trad.mat.dist <- vegdist(trad.mat,method="euclidean")
trad.mat.pcoa <- cmdscale(trad.mat.dist)

# convert and format (pcoa)
trad.pc1.tib <- tibble(label=rownames(trad.mat.pcoa), PC1=trad.mat.pcoa[,1]) %>% 
    mutate(month=str_split_fixed(label,"\\.",6)[,3]) %>%
    mutate(month=fct_relevel(month,c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))) 

# plot
p <- trad.pc1.tib %>%
    ggplot(aes(x=month,y=PC1)) +
        geom_point(color="gray30",size=2.5,alpha=0.8,shape=16) + 
        geom_smooth(aes(x=as.numeric(month)),method="gam",formula=y~s(x,bs="cs",k=4),color="#737145",alpha=0.5,fill="gray90") +
        ggthemes::theme_clean(base_size=12) +
        labs(x="Date", y="Demersal trawl PCoA-1 scores")
#plot(p)
ggsave(filename=here("temp/results/figures/trawl-pcoa.svg"),plot=p,width=130,height=80,units="mm")

# SUPPORTING - plot trawl PCA VARIABLES
p <- plot.PCA(PCA(trad.mat,graph=FALSE,scale.unit=FALSE),choix="var",select="contrib 12",ggoptions=list(size=2,circle.lwd=0.3,circle.lty=3,line.lty=3,line.lwd=0.3)) + theme(plot.title=element_blank(),axis.title=element_text(size=8))
#plot(p)
ggsave(filename=here("temp/results/figures/trawl-pcoa-plot.svg"),plot=p,width=130,height=80,units="mm")
