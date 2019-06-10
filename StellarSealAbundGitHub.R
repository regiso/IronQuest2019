# For IronQuest May 2019
# Steller Sea Lion Abundance  Data


library(agTrend)
library(openxlsx)
library(tidyverse)

# get the raw data from NOAA Fisheries. Look for these files
# file <- "seal.tar.gz"
untar(file, list=TRUE)
untar(file, files="0128190/2.4/data/0-data/SSLnonpupcounts2016.csv")
df1 <- read.csv("0128190/2.4/data/0-data/SSLnonpupcounts2016.csv")


# the lat long for the total 417 sites is at NOAA Fisheres too.
#locs <- "locationslatlong.tar.gz"
untar(locs, list=TRUE)
untar(locs, files="0129877/2.3/data/0-data/NCEI-0129877_US_SSL_Sites.csv")
locs1 <- read.csv("0129877/2.3/data/0-data/NCEI-0129877_US_SSL_Sites.csv")
locsdest <- "C:/Users/OConnor Analytics/Working Files/ironviz/2019/seal/latlong.csv"
locations <- write.csv(locs1, locsdest)

# filter to the breeding months - June & July
# NOAA Fisheries uses these for trend analysis
# Filter data older than 1990 because surveys before then were sparse


# First, use the 2012 data included with the agTrend package.
# The ultimate goal is to duplicate these results with 2016 or 2018 data

data("wdpsNonpups")
wdpsNonpups <- wdpsNonpups %>% filter(YEAR >= 1990) %>%
  droplevels()
wdpsNonpups %>% group_by(SITE) %>% 
  summarise(num_counts=sum(COUNT>0)) %>% ungroup() -> nz_counts
wdpsNonpups %>% left_join(nz_counts) %>% filter(num_counts>1) %>% select(-num_counts) %>% 
  mutate(SITE=factor(SITE)) -> wdpsNonpups
#> Joining, by = "SITE"
wdpsNonpups %>% mutate(obl = as.integer(YEAR<2004)) %>% as.data.frame() -> wdpsNonpups



# The next step is to create prediction and availability models for each site 
# based on the number of surveys.


# The trend (prediction) models will be:
#  0-5  nonzero counts = constant trend signal
#  6-10 nonzero counts = linear trend signal
#  10+  nonzero counts = RW2 (i.e., spline) trend signal)

# The zero inflation (availability) models are:
#  0-5 surveys = constant inflation effect
#  5   surveys = linear inflation effect
#  All surveys have nonzero counts = no availability model




wdpsModels <- wdpsNonpups %>%
  group_by(SITE) %>%
  summarize(
    num_surv = n(),
    nz_count = sum(COUNT>0)) %>%
  ungroup() %>%
  mutate(
    trend = as.character(cut(nz_count, c(0,5,10,30), labels = c("const", "lin", "RW"))),
    avail = as.character(cut(num_surv, c(0,5,30),    labels = c("const", "lin")))
  ) %>%
  mutate(avail = if_else(num_surv==nz_count, "none", avail)) %>%
  select(-num_surv, -nz_count) %>%
  as.data.frame()



# This next code is copied from github and creates prior distributions for the Markov Chain
# Monte Carlo site updating.
# The photo correction adjusts for the count differences observed by the 2 types of 
# photos - overhead and oblique

data("photoCorrection")
photoCorrection <- photoCorrection %>%
  mutate(log_ratio = OBLIQUE/VERTICAL)
gamma.0  <- photoCorrection %>% summarize(mean(log_ratio)) %>% as.double()
gamma.se <- photoCorrection %>% summarize(sd(log_ratio)/sqrt(n())) %>% as.double()

prior.list <- defaultPriorList(trend.limit = 0.2, model.data = wdpsModels,
                               gamma.mean = gamma.0, gamma.prec = 1/gamma.se^2)


# Before running the MCMC we set the upper bounds for predictive counts. 
# We use 3x the maximum count observed at the site.

upper = wdpsNonpups %>% group_by(SITE) %>% summarize(upper=3*max(COUNT)) %>% ungroup()


set.seed(123) 
fit <- mcmc.aggregate(start=1990, end=2030, data=wdpsNonpups, obs.formula=~obl-1, model.data=wdpsModels, 
                      rw.order=list(omega=2), aggregation="REGION",
                      abund.name="COUNT", time.name="YEAR", site.name="SITE", forecast = TRUE,
                      burn=1000, iter=5000, thin=5, prior.list=prior.list, upper=upper, 
                      keep.site.param=TRUE, keep.site.abund=TRUE, keep.obs.param=TRUE)

# Here are the results
fitdat <- fit$aggregation.pred.summary

# Compute trends for 2008 to 2019 and 2029
trend2019 = updateTrend(fit, 2008, 2019, "pred")
trend2029 = updateTrend(fit, 2008, 2029, "pred")

# Change to % growth

growth2019 = mcmc(100*(exp(trend2019[,7:12])-1))
growth2029 = mcmc(100*(exp(trend2029[,7:12])-1))

chart2019 <- data.frame(
  post.median=round(apply(growth2019, 2, median),2),
  HPD.90=round(HPDinterval(growth2019, 0.95),2)
)
chart2029 <- data.frame(
  post.median=round(apply(growth2029, 2, median),2),
  HPD.90=round(HPDinterval(growth2029, 0.95),2)
)


chart2019$region <- rownames(chart2019)
chart2029$region <- rownames(chart2029)


# Make a SITe ID & Join the site ID to the wdpsNonpups

wdpsNonpups2 <- merge(wdpsNonpups, locs1, by.x="SITE", by.y="SITENAME", all.x = TRUE,
                      all.y = TRUE)


dest  <- "C:/Users/OConnor Analytics/Working Files/ironviz/2019/seal/fitdat.csv"
dest2 <- "C:/Users/OConnor Analytics/Working Files/ironviz/2019/seal/wdpsNonpups.csv"
dest3 <- "C:/Users/OConnor Analytics/Working Files/ironviz/2019/seal/wdpsNonpupsMerg.csv"


write_excel_csv(fitdat, dest)
write_excel_csv(wdpsNonpups, dest2)


# link regions to sitename for mapping


# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
