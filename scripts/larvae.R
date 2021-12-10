#!/usr/bin/env Rscript

# load libs and funs
source(here::here("scripts/load-data.R"))

# collapse
trad.collapsed <- trad.master %>% collapse_taxonomy(rmfw=TRUE,lifestage="Larvae",collapse=TRUE)# LARVAE
edna.collapsed <- edna.filt %>% collapse_taxonomy(rmfw=TRUE,lifestage="Larvae",collapse=TRUE)# LARVAE

# get sample n
glue("\nNumber eDNA samples total = {edna.collapsed %>% filter(partnerID=='MBA' & !grepl('WHIT',eventID)) %>% distinct(sampleHash) %>% count() %>% pull()}",.trim=FALSE)

# get time duration
dur.edna <- edna.collapsed %>% filter(partnerID=="MBA" & !grepl("WHIT",eventID)) %>% distinct(eventDate) %>% arrange(eventDate) %>% pull(eventDate)
# get n samples during duration
glue("\nNumber unique ichthyoplankton samples concurrent with eDNA survey = {
    trad.collapsed %>% 
    filter(partnerID=='MBA' & !grepl('WHIT',eventID)) %>%
    filter(eventDate %within% lubridate::interval(start=first(dur.edna),end=last(dur.edna))) %>%
    distinct(eventID,eventDate,fieldNumber) %>%
    count() %>% 
    pull()
    }",.trim=FALSE)


# make a list of species common to edna and trad (by partnerID)
comb.spp <- full_join(distinct(edna.collapsed,partnerID,species),distinct(trad.collapsed,partnerID,species), by=c("species","partnerID")) %>% arrange(partnerID,species)

# expand for all spp and get groups
trad.expanded <- trad.collapsed %>% expand_and_summarise(sppdf=comb.spp,eventdf=events,primerdf=NA,method="traditional",spatialgroup="partner",temporalgroup="month",correct=NA)# LARVAE
edna.expanded <- edna.collapsed %>% expand_and_summarise(sppdf=comb.spp,eventdf=events,primerdf=primer.bias,method="edna",spatialgroup="partner",temporalgroup="month",correct=FALSE)# LARVAE

# join and clean 
surveys.joined <- join_and_clean(trad=trad.expanded,edna=edna.expanded)

# filter the required data and change the month label
surveys.joined %<>% mutate(month=lubridate::month(temporalGroup,label=TRUE)) %>% filter(partnerID=="MBA" & localityID!="WHIT") # LARVAE/EGGS/SPAWNING

# first run rank and zeroes filters to see data and set cutoff
glue("\nSpecies ranked by total indiv abundance ...",.trim=FALSE)
surveys.joined %>% rank_and_subset(rankmethod="individuals",print=TRUE,n="all") %>% invisible() # LARVAE
glue("\nSpecies ranked by proportion zero catch events ...",.trim=FALSE)
surveys.joined %>% zero_catch_filter(comparison="individuals",filter="both0",print=TRUE,maxprop0=1) %>% invisible()

# subset by ranked indivs 
#surveys.joined %<>% rank_and_subset(rankmethod="individuals",print=FALSE,n=26)# LARVAE

# filter rare species with few records
surveys.joined %<>% zero_catch_filter(comparison="individuals",filter="both0",print=FALSE,maxprop0=0.99)

# get number samples
glue("\nData included in model ...",.trim=FALSE)
glue("Total number species = {surveys.joined %>% distinct(species) %>% count() %>% pull(n)}")
glue("Total number samples = {surveys.joined %>% distinct(temporalGroup,localitySite,sampleHash) %>% count() %>% pull(n)}")
glue("Total number sampling dates = {surveys.joined %>% distinct(temporalGroup) %>% count() %>% pull(n)}")

# reduce LARVAE
surveys.joined.coll <- surveys.joined %>% 
    distinct(species,temporalGroup,spatialGroup,readsByGroupProportion,individualsByGroupRate)

# summarise
surveys.joined.avg <- surveys.joined.coll %>% 
    group_by(species) %>%
    summarise(readsByGroupProportionMean=mean(fourth_root(readsByGroupProportion),na.rm=TRUE), sem=se(fourth_root(readsByGroupProportion)), individualsByGroupRateMean=mean(fourth_root(individualsByGroupRate),na.rm=TRUE), .groups="drop")

# check lms
pdf(file=here("temp/results/figures/larvae-check-lm.pdf"))
    lm(readsByGroupProportionMean ~ individualsByGroupRateMean, data=surveys.joined.avg) %>% performance::check_model() %>% plot()
dev.off()

glue("\nLinear model summary ...",.trim=FALSE)
lm(readsByGroupProportionMean ~ individualsByGroupRateMean, data=surveys.joined.avg) %>% summary()

# plot LARVAE
p <- surveys.joined.avg %>% ggplot(aes(y=readsByGroupProportionMean,x=individualsByGroupRateMean,ymin=readsByGroupProportionMean-sem,ymax=readsByGroupProportionMean+sem)) + 
        geom_pointrange(color="gray30",size=0.5) +
        geom_text(aes(label=species),check_overlap=TRUE,vjust=1.5,size=3,color="gray50") +
        annotate(geom="label",x=1.25,y=0.05,label=extract_p(surveys.joined.avg,y="readsByGroupProportionMean",x="individualsByGroupRateMean",type="lm",dp=9),size=3) +
        geom_smooth(method="lm",formula=y~x,alpha=0.5,color="#2f8685",fill="gray90") +
        theme_clean(base_size=12) +
        labs(x="Larval ichthyoplankton abundance\n(4th root transformed CPUE)",y="Proportion of fish community (eDNA)\n(4th root  transformed CPUE)")
#plot(p)
ggsave(filename=here("temp/results/figures/larvae.svg"),plot=p,width=120,height=120,units="mm")

# report
glue("\nFigures saved to 'temp/results/figures'",.trim=FALSE)
