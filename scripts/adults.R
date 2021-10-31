#!/usr/bin/env Rscript

# load libs and funs
source(here::here("scripts/load-data.R"))

# collapse taxonomy
trad.collapsed <- trad.master %>% collapse_taxonomy(rmfw=TRUE,lifestage="Adult",collapse=TRUE)# ADULTS/SPAWNING
edna.collapsed <- edna.filt %>% collapse_taxonomy(rmfw=TRUE,lifestage="Adult",collapse=TRUE)# ADULTS/SPAWNING

# get sample n
glue("\nNumber L4 eDNA samples total = {edna.collapsed %>% filter(partnerID=='MBA' & grepl('L4',eventID)) %>% distinct(sampleHash) %>% count() %>% pull()}",.trim=FALSE)

# get time duration
dur.edna <- edna.collapsed %>% filter(partnerID=="MBA" & grepl("L4",eventID)) %>% distinct(eventDate) %>% arrange(eventDate) %>% pull(eventDate)
# get n samples during duration
glue("\nNumber unique L4 trawls concurrent with eDNA survey = {
    trad.collapsed %>% 
    filter(partnerID=='MBA' & grepl('L4',eventID)) %>%
    filter(eventDate %within% lubridate::interval(start=first(dur.edna),end=last(dur.edna))) %>%
    distinct(eventID,eventDate,fieldNumber) %>%
    count() %>% 
    pull()
    }",.trim=FALSE)

# make a list of species common to edna and trad (by partnerID)
comb.spp <- full_join(distinct(edna.collapsed,partnerID,species),distinct(trad.collapsed,partnerID,species),by=c("species","partnerID")) %>% arrange(partnerID,species)

# expand the tables
trad.expanded <- trad.collapsed %>% expand_and_summarise(sppdf=comb.spp,eventdf=events,primerdf=NA,method="traditional",spatialgroup="locality",temporalgroup="day",correct=NA)# EGGS/ADULTS/SPAWNING
edna.expanded <- edna.collapsed %>% expand_and_summarise(sppdf=comb.spp,eventdf=events,primerdf=primer.bias,method="edna",spatialgroup="locality",temporalgroup="day",correct=FALSE)# EGGS/ADULTS/SPAWNING

# join and clean 
surveys.joined <- join_and_clean(trad=trad.expanded,edna=edna.expanded)

# filter the required data and change the month label
surveys.joined %<>% mutate(month=lubridate::month(temporalGroup,label=TRUE)) %>% filter(partnerID=="MBA" & localityID=="L4") # ADULTS

# first run rank and zeroes filters to see data and set cutoff
glue("\nSpecies ranked by total eDNA abundance ...",.trim=FALSE)
surveys.joined %>% rank_and_subset(rankmethod="reads",print=TRUE,n="all") %>% invisible()
glue("\nSpecies ranked by proportion zero catch events ...",.trim=FALSE)
surveys.joined %>% zero_catch_filter(comparison="individuals",filter="both0",print=TRUE,maxprop0=1) %>% invisible()

# subset by ranked indivs 
surveys.joined %<>% rank_and_subset(rankmethod="reads",print=FALSE,n=63)# @ 1000 reads

# filter rare species with few records
surveys.joined %<>% zero_catch_filter(comparison="individuals",filter="both0",print=FALSE,maxprop0=0.90)# previously 0.95

# load up the fishbase data and merge with joined surveys 
tax.table.ann <- tax.table %>% annotate_fishbase(collapse="Adult")# ADULTS

# join
surveys.joined.ann <- surveys.joined %>% left_join(tax.table.ann,by=c("month","species"))

# print and check the data
# generate some stats
glue("\neDNA data included in model ...",.trim=FALSE)
glue("Total number species = {surveys.joined.ann %>% distinct(species) %>% count() %>% pull(n)}")
glue("Total number samples = {surveys.joined.ann %>% distinct(sampleHash) %>% count() %>% pull(n)}")
glue("Total number sampling dates = {surveys.joined.ann %>% distinct(temporalGroup) %>% count() %>% pull(n)}")

# collapse
surveys.joined.coll <- surveys.joined.ann %>% 
    distinct(species,lifestyle,temporalGroup,readsByGroupProportion,individualsByGroupRate)

# summarise and transform
surveys.joined.avg <- surveys.joined.coll %>% 
    group_by(species) %>%
    summarise(readsByGroupProportionMean=mean(fourth_root(readsByGroupProportion),na.rm=TRUE), sem=se(fourth_root(readsByGroupProportion)), individualsByGroupRateMean=mean(fourth_root(individualsByGroupRate),na.rm=TRUE), .groups="drop")

# test lms
pdf(file=here("temp/results/figures/adults-check-lm.pdf"))
    lm(readsByGroupProportionMean ~ individualsByGroupRateMean, data=surveys.joined.avg) %>% performance::check_model() %>% plot()
dev.off()

glue("\nLinear model summary ...",.trim=FALSE)
lm(readsByGroupProportionMean ~ individualsByGroupRateMean, data=surveys.joined.avg) %>% summary()

# plot ADULTS
p <- surveys.joined.avg %>% ggplot(aes(y=readsByGroupProportionMean,x=individualsByGroupRateMean,ymin=readsByGroupProportionMean-sem,ymax=readsByGroupProportionMean+sem)) + 
        geom_pointrange(color="gray30",size=0.5) +
        geom_text(aes(label=species),check_overlap=TRUE,vjust=1.5,size=3,color="gray50") +
        annotate(geom="label",x=5,y=0.05,label=extract_p(df=surveys.joined.avg,y="readsByGroupProportionMean",x="individualsByGroupRateMean",type="lm",dp=7),size=3) +
        geom_smooth(method="lm",formula=y~x,alpha=0.5,color="#737145",fill="gray90") +
        theme_clean(base_size=12) +
        labs(x="Demersal trawl adult abundance\n(4th root transformed CPUE)",y="Proportion of fish community (eDNA)\n(4th root transformed CPUE)")
#plot(p)
ggsave(filename=here("temp/results/figures/adults.svg"),plot=p,width=120,height=120,units="mm")


# plot RESIDUALS
p <- surveys.joined.coll %>% 
    lm(formula=fourth_root(readsByGroupProportion)~fourth_root(individualsByGroupRate),data=.,na.action=na.omit) %>% #summary
    broom::augment(data=surveys.joined.coll,type.predict="response") %>%
    group_by(lifestyle) %>%
    summarise(mean.resid=mean(.resid,na.rm=TRUE),se.resid=se(.resid),.groups="drop") %>%
    mutate(lifestyle=fct_relevel(lifestyle,"pelagic","benthopelagic","benthic")) %>%
    ggplot(aes(y=mean.resid,x=fct_reorder(lifestyle,mean.resid),ymin=mean.resid-se.resid,ymax=mean.resid+se.resid,color=lifestyle)) +
        geom_pointrange(size=1) + 
        geom_text(aes(label=lifestyle),vjust=-2,color="gray50") +
        geom_hline(yintercept=0,linetype=2) +
        theme_clean(base_size=12) +
        scale_color_ptol(name=NULL) + 
        theme(axis.text.x=element_blank(),legend.position="none") +
        labs(x="Life history",y="Index of eDNA abundance relative to trawl catch\n(linear model residual values)")
#plot(p)
ggsave(filename=here("temp/results/figures/residuals-lm.svg"),plot=p,width=90,height=120,units="mm")


# run glmmTMB - NEGATIVE BINOMIAL MODEL

# subset dataframe and reformat
surveys.joined.ann.sub <- surveys.joined.ann %>% 
    select(nReads,readsBySampleTotal,individualsByGroupRate,maxEfficiency,spawningByMonthUK,lifestyle,temporalGroup,localitySite,sampleHash,order,family,species) %>% 
    mutate_if(~!is.numeric(.),factor)

# report
glue("\nRunning negative binomial model - may take a few minutes ...",.trim=FALSE)

# run model
m0 <- glmmTMB(
    formula=nReads~offset(log(readsBySampleTotal))+scale(individualsByGroupRate)+scale(maxEfficiency)+spawningByMonthUK+lifestyle+(1|temporalGroup/localitySite/sampleHash)+(1|order/family/species),
    data=surveys.joined.ann.sub, 
    family=nbinom2(link="log"), 
    ziformula=~.,
    REML=FALSE, 
    na.action=na.omit)

# drop the spawning month
m1 <- update(m0,.~. -spawningByMonthUK)


# summary
glue("\nPrinting model summary in standard format ...\n",.trim=FALSE)
summary(m0)
#m0 %>% broom.mixed::tidy()

# tabulate model output
glue("\nPrinting model summary in LaTeX format ...\n",.trim=FALSE)
m0 %>% 
    broom.mixed::tidy() %>% 
    mutate_if(is.character,~str_replace_all(.,"_","")) %>%
    xtable(caption="blahhh",display=c("s","s","s","s","s","f","f","f","e")) %>% 
    print.xtable(include.rownames=FALSE,booktabs=TRUE,sanitize.text.function=identity,caption.placement="top",size="scriptsize")

# model checks
glue("\nRunning model checks - may take a few minutes ...",.trim=FALSE)
 
# run DHARMa residual sims on m0
sim <- simulateResiduals(fittedModel=m0, plot=FALSE, n=1000, seed=42)

pdf(file=here("temp/results/figures/model-check-dharma.pdf"))
    # individual fit tests
    suppressMessages({
        testUniformity(sim,plot=TRUE)# tests if the overall distribution conforms to expectations
        testOutliers(sim,plot=TRUE)# tests if there are more simulation outliers than expected
        testDispersion(sim,plot=TRUE)# tests if the simulated dispersion is equal to the observed dispersion
        testQuantiles(sim,plot=TRUE)# fits a quantile regression or residuals against a predictor (default predicted value), and tests of this conforms to the expected quantile
        testZeroInflation(sim,plot=TRUE)# tests if there are more zeros than expected
    })
dev.off()

# run performance checks
perm <- performance::check_model(m0)

# plot pdf
pdf(file=here("temp/results/figures/model-check-performance.pdf"),height=12,width=18)
    plot(perm)
dev.off()

# model checks
glue("\nModel selection and estimated marginal means ...\n",.trim=FALSE)

# compare performance of models
#performance::compare_performance(m0,m1,rank=TRUE)
performance::r2_zeroinflated(m0,method="correlation")
performance::r2_zeroinflated(m1,method="correlation")
glue("\n",.trim=FALSE)
MuMIn::model.sel(m0,m1)

# get marginal means and ratios
#emmeans(m0,"spawningByMonthUK",type="response")
#emmeans(m0,"spawningByMonthUK",type="response") %>% plot()
glue("\n",.trim=FALSE)
pairs(emmeans(m0,"spawningByMonthUK"),type="response",reverse=TRUE)

# report
glue("\nFigures saved to 'temp/results/figures'",.trim=FALSE)
