#!/usr/bin/env Rscript

# load libs and funs
source(here::here("scripts/load-data.R"))

# collapse
trad.collapsed <- trad.master %>% collapse_taxonomy(rmfw=TRUE,lifestage="Adult",collapse=TRUE)# ADULTS/SPAWNING
edna.collapsed <- edna.filt %>% collapse_taxonomy(rmfw=TRUE,lifestage="Adult",collapse=TRUE)# ADULTS/SPAWNING

# make a list of species common to edna and trad (by partnerID)
comb.spp <- full_join(distinct(edna.collapsed,partnerID,species),distinct(trad.collapsed,partnerID,species)) %>% arrange(partnerID,species)

# expand
trad.expanded <- trad.collapsed %>% expand_and_summarise(sppdf=comb.spp,eventdf=events,primerdf=NA,method="traditional",spatialgroup="locality",temporalgroup="day",correct=NA)# EGGS/ADULTS/SPAWNING
edna.expanded <- edna.collapsed %>% expand_and_summarise(sppdf=comb.spp,eventdf=events,primerdf=primer.bias,method="edna",spatialgroup="locality",temporalgroup="day",correct=FALSE)# EGGS/ADULTS/SPAWNING

# join and clean 
surveys.joined <- join_and_clean(trad=trad.expanded,edna=edna.expanded)

# filter the required data and change the month label
surveys.joined %<>% mutate(month=lubridate::month(temporalGroup,label=TRUE)) %>% filter(partnerID=="MBA" & localityID!="WHIT") # LARVAE/EGGS/SPAWNING

# first run rank and zeroes filters to see data and set cutoff
surveys.joined %>% rank_and_subset(rankmethod="reads",print=TRUE,n="all") %>% invisible() # ADULTS/SPAWNING
surveys.joined %>% zero_catch_filter(comparison="individuals",filter="both0",print=TRUE,maxprop0=1) %>% invisible()

# subset by ranked indivs 
surveys.joined %<>% rank_and_subset(rankmethod="reads",print=FALSE,n=75)

# filter rare species with few records
surveys.joined %<>% zero_catch_filter(comparison="individuals",filter="both0",print=FALSE,maxprop0=0.95)

# load up the fishbase data and merge with joined surveys 
tax.table.ann <- tax.table %>% annotate_fishbase(collapse="Adult")# ADULTS

# join
surveys.joined.ann <- surveys.joined %>% left_join(tax.table.ann)

# print and check the data
print(surveys.joined.ann,n=20)
surveys.joined.ann %>% distinct(species) %>% print(n=Inf)
surveys.joined.ann %>% distinct(sampleHash) %>% print(n=Inf)

# reduce 
surveys.joined.coll <- surveys.joined.ann %>% mutate(binvar=if_else(spawningByMonthUK==TRUE,1,0))

# plot SPAWNING GLMs
p <-surveys.joined.coll %>% ggplot(aes(x=nReads,y=binvar)) +
            geom_smooth(method="glm",method.args=list(family=binomial(link="logit"),na.action=na.omit),se=TRUE,alpha=0.5,color="#737145",fill="gray90") +
            geom_rug(data=filter(surveys.joined.coll,binvar==0),sides="b",alpha=0.2) +
            geom_rug(data=filter(surveys.joined.coll,binvar==1),sides="t",alpha=0.2) +
            annotate(geom="label",x=90000,y=0.1,label=extract_p(surveys.joined.coll,y="binvar",x="nReads",type="glm",dp=16),size=3) +
            labs(x="eDNA read abundance",y="Breeding by month (Pr)") +
            ggthemes::theme_clean(base_size=12)
plot(p)
ggsave(filename=here("temp/results/figures/spawning.svg"),plot=p,width=110,height=140,units="mm")
