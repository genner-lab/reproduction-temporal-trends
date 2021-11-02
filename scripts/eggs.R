#!/usr/bin/env Rscript

# load libs and funs
source(here::here("scripts/load-data.R"))

# collapse eggs
trad.collapsed <- trad.master %>% collapse_taxonomy(rmfw=TRUE,lifestage="Eggs",collapse=TRUE)# EGGS
edna.collapsed <- edna.filt %>% collapse_taxonomy(rmfw=TRUE,lifestage="Eggs",collapse=TRUE)# EGGS

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
comb.spp <- full_join(distinct(edna.collapsed,partnerID,species),distinct(trad.collapsed,partnerID,species),by=c("species","partnerID")) %>% arrange(partnerID,species)

# expand
trad.expanded <- trad.collapsed %>% expand_and_summarise(sppdf=comb.spp,eventdf=events,primerdf=NA,method="traditional",spatialgroup="locality",temporalgroup="day",correct=NA)# EGGS
edna.expanded <- edna.collapsed %>% expand_and_summarise(sppdf=comb.spp,eventdf=events,primerdf=primer.bias,method="edna",spatialgroup="locality",temporalgroup="day",correct=FALSE)# EGGS

# join and clean 
surveys.joined <- join_and_clean(trad=trad.expanded,edna=edna.expanded)

# filter the required data and change the month label
surveys.joined %<>% mutate(month=lubridate::month(temporalGroup,label=TRUE)) %>% filter(partnerID=="MBA" & localityID!="WHIT") # EGGS

# get stats
glue("\nData included in model ...",.trim=FALSE)
glue("Total number species = {surveys.joined %>% distinct(species) %>% count() %>% pull(n)}")
glue("Total number samples = {surveys.joined %>% distinct(temporalGroup,localitySite,sampleHash) %>% count() %>% pull(n)}")
glue("Total number sampling dates = {surveys.joined %>% distinct(temporalGroup) %>% count() %>% pull(n)}")

# reduce and plot EGGS
surveys.joined.coll <- surveys.joined %>% 
    select(species,sampleHash,temporalGroup,spatialGroup,readsBySampleProportion,individualsByGroupRate) %>%
    filter(species=="Sardina pilchardus") %>% 
    group_by(temporalGroup,spatialGroup) %>%
    summarise(readsByGroupProportionMean=mean(fourth_root(readsBySampleProportion),na.rm=TRUE), sem=se(fourth_root(readsBySampleProportion)), individualsByGroupRateMean=mean(individualsByGroupRate,na.rm=TRUE),.groups="drop")

# check lms
# check lms
pdf(file=here("temp/results/figures/larvae-check-lm.pdf"))
    lm(readsByGroupProportionMean ~ individualsByGroupRateMean, data=surveys.joined.coll) %>% performance::check_model() %>% plot()
dev.off()

glue("\nLinear model summary ...",.trim=FALSE)
lm(readsByGroupProportionMean ~ individualsByGroupRateMean, data=surveys.joined.coll) %>% summary()

# plot EGGS
p <- surveys.joined.coll %>% ggplot(aes(y=readsByGroupProportionMean,x=individualsByGroupRateMean,ymin=readsByGroupProportionMean-sem,ymax=readsByGroupProportionMean+sem)) + 
        geom_pointrange(color="gray30",size=0.5) +
        annotate(geom="label",x=225,y=0.1,label=extract_p(surveys.joined.coll,y="readsByGroupProportionMean",x="individualsByGroupRateMean",type="lm",dp=3),size=3) +
        geom_smooth(method="lm",formula=y~x,alpha=0.5,color="#2f8685",fill="gray90") +
        theme_clean(base_size=12) +
        labs(x="Pilchard egg ichthyoplankton abundance\n(CPUE)",y="Proportion of fish community (eDNA)\n(4th root transformed CPUE")
#plot(p)
ggsave(filename=here("temp/results/figures/eggs.svg"),plot=p,width=120,height=120,units="mm")

# report
glue("\nFigures saved to 'temp/results/figures'",.trim=FALSE)
