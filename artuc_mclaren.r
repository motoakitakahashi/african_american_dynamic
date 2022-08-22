# Motoaki Takahashi
# May 2022
# R 4.2.0

rm(list=ls())

# for the main desktop
setwd("D:/onedrive/OneDrive - The Pennsylvania State University/dynamic/mig_wage")

# for my laptop
# setwd("C:/Users/takah/OneDrive - The Pennsylvania State University/dynamic/mig_wage")

library(FENmlm)
library(fixest)
library(haven)

library(tidyr)
library(dplyr)
library(readxl)
library(readr)
library(data.table)
library(haven)
library(stargazer)
library(ggplot2)
library(jtools)
library(ggstance)
library(broom.mixed)
library(stringr)
library(openxlsx)
library(texreg)
library(geosphere)
library(sp)
library(survey)
library(ipumsr)
library(primes)
library(ivreg)

# read data
load("output/r/migwagerent.rda")



# set parameter values
N = length(unique(data$my_geo)) # number of locations

C = 8 # the periods of life

gamma1 = 0.25 # Cobb-douglas share on housing



geolist = sort(unique(data$my_geo))

# assign prime numbers to my_geo
geopri <- generate_n_primes(N)

geolist_df = data.frame(geolist, geopri)

names(geolist_df) <- c("my_geo", "my_geopri")



data = left_join(data, geolist_df, by = "my_geo")


origin_df = geolist_df
names(origin_df) <- c("origin", "originpri")

data = left_join(data, origin_df, by = "origin")

# product of originpri and my_geopri uniquely identify an unordered pair of locations

data <- mutate(data, unop = originpri * my_geopri)


data_loc = select(data, race2, year, age, coh, my_geo)
data_loc = unique(data_loc)

# make a similar object for regression without staying
data_loc2 = data_loc

data = mutate(data, twodec = floor((year-1940)/20)) # two decades
M = max(data$unop)
data = mutate(data, own = (M+1)*(origin == my_geo))
data = mutate(data, cons = 1)
data = mutate(data, unop2 = unop * (unop > own))
# unop2 signifies all unordered pairs of "different" locations
# unop2 sets 0 for a pair of the same locations

####################################################################################

# the first step without staying (migrating from a location to the same location)

# drop observations for staying
data2  = data %>%
  mutate(own = my_geo / origin) %>%
  filter(own != 1)

data2 = select(data2, -own)

data2 = mutate(data2, cons = 1)

# make dummy variables for South-to-North migration

# read the concordance table that maps my geographic units to the North/South
geosouth = read_csv("../geo_unit/output/my_geo_south.csv")
names(geosouth)[2] = "dessouth"

orisouth = geosouth
names(orisouth) = c("origin", "orisouth")

data2 = left_join(data2, geosouth, by = "my_geo")
data2 = left_join(data2, orisouth, by = "origin")

data2 = data2 %>%
  mutate(oridess = paste(orisouth, dessouth, sep = ""))

##################################################################################
# when I need to rerun this code, rerun the entire code except the following
# first step without staying
first_step_oridess = fepois(mig ~ 1 | cons + age + year + race2 + year^age^race2^origin +
                       year^age^race2^my_geo + unop + race2 ^ unop + twodec ^ unop +
                      age ^ unop + race2 ^ twodec ^ unop +oridess,
                    data = data2, combine.quick=FALSE)
# race2 ^ twodec ^ unop + 
save(first_step_oridess, file="output/r/mig_elas_first_step_oridess.rda")
###################################################################################
load(file="output/r/mig_elas_first_step_oridess.rda")


summary(first_step_oridess)
first_step_oridess_fe = fixef(first_step_oridess)
first_step_oridess_coef = coef(first_step_oridess)

first_step_oridess_p = predict(first_step_oridess, data2)

# dataframe to store migration costs (or estimated fixed effects in the first step)
mig_cost_df2 = select(data2, race2, year, age, my_geo, origin, my_geopri, originpri, unop, twodec, cons, oridess)


# cons
cons_fe2 = data.frame(as.integer(names(first_step_oridess_fe$cons)), first_step_oridess_fe$cons)
names(cons_fe2) = c("cons", "cons_fe")

# age
age_fe2 = data.frame(as.integer(names(first_step_oridess_fe$age)), first_step_oridess_fe$age)
names(age_fe2) = c("age", "age_fe")

# year
year_fe2 = data.frame(as.integer(names(first_step_oridess_fe$year)), first_step_oridess_fe$year)
names(year_fe2) = c("year", "year_fe")

# race2
race2_fe2 = data.frame(as.integer(names(first_step_oridess_fe$race2)), first_step_oridess_fe$race2)
names(race2_fe2) = c("race2", "race2_fe")

# unordered pairs of locations (unop)
unop_fe2 = data.frame(as.integer(names(first_step_oridess_fe$unop)), first_step_oridess_fe$unop)
names(unop_fe2) = c("unop", "unop_fe")

# oridess
oridess_fe2 = data.frame(as.integer(names(first_step_oridess_fe$oridess)), first_step_oridess_fe$oridess)
names(oridess_fe2) = c("oridess", "oridess_fe")

# race2 and unop
race2_unop_fe2 = data.frame(names(first_step_oridess_fe$`race2^unop`), first_step_oridess_fe$`race2^unop`)
names(race2_unop_fe2) = c("race2_unop", "race2_unop_fe")
split_race2_unop2 = str_split(race2_unop_fe2$race2_unop , "_", simplify = TRUE)

split_race2_unop2 = as.data.frame(split_race2_unop2)
names(split_race2_unop2) = c("race2", "unop")
split_race2_unop2$race2 = as.integer(split_race2_unop2$race2)
split_race2_unop2$unop = as.integer(split_race2_unop2$unop)

race2_unop_fe2 = cbind(race2_unop_fe2, split_race2_unop2)
race2_unop_fe2 = select(race2_unop_fe2, -race2_unop)

rm(split_race2_unop2)

# two decades (twodec) and unop
twodec_unop_fe2 = data.frame(names(first_step_oridess_fe$`twodec^unop`), first_step_oridess_fe$`twodec^unop`)
names(twodec_unop_fe2) = c("twodec_unop", "twodec_unop_fe")
split_twodec_unop2 = str_split(twodec_unop_fe2$twodec_unop , "_", simplify = TRUE)
split_twodec_unop2 = as.data.frame(split_twodec_unop2)

names(split_twodec_unop2) = c("twodec", "unop")
split_twodec_unop2$twodec = as.integer(split_twodec_unop2$twodec)
split_twodec_unop2$unop = as.integer(split_twodec_unop2$unop)

twodec_unop_fe2 = cbind(twodec_unop_fe2, split_twodec_unop2)
twodec_unop_fe2 = select(twodec_unop_fe2, -twodec_unop)

rm(split_twodec_unop2)

# race2-twodec-unop
race2_twodec_unop_fe2 = data.frame(names(first_step_oridess_fe$`race2^twodec^unop`), first_step_oridess_fe$`race2^twodec^unop`)
names(race2_twodec_unop_fe2) = c("race2_twodec_unop", "race2_twodec_unop_fe")
split_race2_twodec_unop2 = str_split(race2_twodec_unop_fe2$race2_twodec_unop , "_", simplify = TRUE)
split_race2_twodec_unop2 = as.data.frame(split_race2_twodec_unop2)

names(split_race2_twodec_unop2) = c("race2", "twodec", "unop")
split_race2_twodec_unop2$race2 = as.integer(split_race2_twodec_unop2$race2)
split_race2_twodec_unop2$twodec = as.integer(split_race2_twodec_unop2$twodec)
split_race2_twodec_unop2$unop = as.integer(split_race2_twodec_unop2$unop)

race2_twodec_unop_fe2 = cbind(race2_twodec_unop_fe2, split_race2_twodec_unop2)
race2_twodec_unop_fe2 = select(race2_twodec_unop_fe2, -race2_twodec_unop)

rm(split_race2_twodec_unop2)

# age-unop
age_unop_fe2 = data.frame(names(first_step_oridess_fe$`age^unop`), first_step_oridess_fe$`age^unop`)
names(age_unop_fe2) = c("age_unop", "age_unop_fe")
split_age_unop2 = str_split(age_unop_fe2$age_unop , "_", simplify = TRUE)
split_age_unop2 = as.data.frame(split_age_unop2)

names(split_age_unop2) = c("age", "unop")
split_age_unop2$age = as.integer(split_age_unop2$age)
split_age_unop2$unop = as.integer(split_age_unop2$unop)

age_unop_fe2 = cbind(age_unop_fe2, split_age_unop2)
age_unop_fe2 = select(age_unop_fe2, -age_unop)

rm(split_age_unop2)

# collect fixed effects in a single dataframe
mig_cost_df2 = left_join(mig_cost_df2, cons_fe2, by = "cons")

mig_cost_df2$oridess = as.integer(mig_cost_df2$oridess)
mig_cost_df2 = left_join(mig_cost_df2, oridess_fe2, by = "oridess")
mig_cost_df2 = left_join(mig_cost_df2, age_fe2, by = "age")
mig_cost_df2 = left_join(mig_cost_df2, year_fe2, by = "year")
mig_cost_df2 = left_join(mig_cost_df2, race2_fe2, by = "race2")
mig_cost_df2 = left_join(mig_cost_df2, unop_fe2, by = "unop")
mig_cost_df2 = left_join(mig_cost_df2, race2_unop_fe2, by = c("race2", "unop"))
mig_cost_df2 = left_join(mig_cost_df2, twodec_unop_fe2, by = c("twodec", "unop"))
mig_cost_df2 = left_join(mig_cost_df2, race2_twodec_unop_fe2, by = c("race2", "twodec", "unop"))
mig_cost_df2 = left_join(mig_cost_df2, age_unop_fe2, by = c("age", "unop"))

mig_cost_df2 = mutate(mig_cost_df2, tau_tilde = cons_fe + 
                        age_fe + year_fe + race2_fe + unop_fe + race2_unop_fe +
                        twodec_unop_fe + race2_twodec_unop_fe + age_unop_fe + oridess_fe)


# year-age-race2-origin
first_step_oridess_fe$`year^age^race2^origin`

yaro_df2 = data.frame(names(first_step_oridess_fe$`year^age^race2^origin`), first_step_oridess_fe$`year^age^race2^origin`)

names(yaro_df2) = c("yaro", "yarofe")

yaro_df2 = yaro_df2 %>%  mutate(year = substr(yaro, 1, 4)) %>% mutate(age = substr(yaro, 6, 7)) %>% mutate(race2 = substr(yaro, 9, 9)) %>% 
  mutate(origin = substring(yaro, 11) )
yaro_df2$year = as.integer(yaro_df2$year)
yaro_df2$race2 = as.integer(yaro_df2$race2)
yaro_df2$origin = as.integer(yaro_df2$origin)
yaro_df2$age = as.integer(yaro_df2$age)

yaro_df2 = mutate(yaro_df2, year = year - 10)
yaro_df2 = mutate(yaro_df2, age = age - 10)

yaro_df2 = select(yaro_df2, -yaro)

names(yaro_df2)[5] = "my_geo"

data_loc2 = left_join(data_loc2, yaro_df2, by = c("my_geo", "year", "race2", "age"))

# year-race2-my_geo
first_step_oridess_fe$`year^age^race2^my_geo`

yard_df2 = data.frame(names(first_step_oridess_fe$`year^age^race2^my_geo`), first_step_oridess_fe$`year^age^race2^my_geo`)
names(yard_df2) = c("yard", "yardfe")

yard_df2 = yard_df2 %>% mutate(year = substr(yard, 1, 4)) %>% mutate(age = substr(yard, 6, 7)) %>% 
  mutate(race2 = substr(yard, 9, 9)) %>% mutate(my_geo = substring(yard, 11) )
yard_df2 = select(yard_df2, -yard)

yard_df2$year = as.integer(yard_df2$year)
yard_df2$race2 = as.integer(yard_df2$race2)
yard_df2$my_geo = as.integer(yard_df2$my_geo)
yard_df2$age = as.integer(yard_df2$age)

data_loc2 = left_join(data_loc2, yard_df2, by = c("year", "my_geo", "race2", "age"))

# make variables for the second step
# temporarily assume that all s are one
data2 = mutate(data2, rw = adj_pyrl_pc_r_c_i_2010/(mrmg2010)^gamma1)
data2 = mutate(data2, logw = log(adj_pyrl_pc_r_c_i_2010), logrw = log(rw))

# make current real wages
wagedf2 = select(data2, year, age, race2, my_geo, logw, logrw)

wagedf2 = unique(wagedf2)

data_loc2 = left_join(data_loc2, wagedf2, by = c("year", "age", "race2", "my_geo"))

# merge data_lo with the survival probability data
load("../life_table/output/sall.rda")

data_loc2 = left_join(data_loc2, sp_all_long, by = c("age", "year", "race2"))

# make lagged real wages
laggedwagedf2 = select(data2, year, age, race2, my_geo, logw, logrw)
laggedwagedf2 = laggedwagedf2 %>% mutate(year = year + 10)
names(laggedwagedf2)[c(5, 6)] = c("laglogw", "laglogrw")
laggedwagedf2 = unique(laggedwagedf2)

data_loc2 = left_join(data_loc2, laggedwagedf2, by = c("year", "age", "race2", "my_geo"))

#################################################################################
# The second step starts here



# (2.2) the second step using the first step without staying

# make the explained variable in the second step
data_loc2 = mutate(data_loc2, lhs = yardfe + s * yarofe)

# make the explanatory variables
data_loc2 = mutate(data_loc2, slogw = s*logw, slogrw = s*logrw)
data_loc2 = mutate(data_loc2, slaglogw = s*laglogw, slaglogrw = s*laglogrw)

# the new explained variable
data_loc2 = mutate(data_loc2, lhs_new = (yardfe/s) + yarofe)

# new regression
fe_iv_rw3_n2 = feols(lhs_new ~ 1|year + race2 + my_geo + age +race2^my_geo + age^my_geo +year^my_geo + 
                      year^race2 + race2^my_geo^year  |
                      logrw ~ laglogrw, data = data_loc2, cluster = c("my_geo"))
summary(fe_iv_rw3_n2)

fe_iv_rw3_n_h2 = feols(lhs_new ~ 1|year + race2 + my_geo + age +race2^my_geo + age^my_geo +year^my_geo + 
                        year^race2 + race2^my_geo^year  |
                        logrw ~ laglogrw, data = data_loc2, vcov = "hetero")
summary(fe_iv_rw3_n_h2)

fe_iv_rw6_n2 = feols(lhs_new ~ 1|year + race2 + my_geo + age +race2^my_geo + age^my_geo +year^my_geo + 
                      age^race2 + year^race2 +race2^my_geo^year + race2^my_geo^age  |
                      logrw ~ laglogrw, data = data_loc2, cluster = c("my_geo"))
summary(fe_iv_rw6_n2)

fe_iv_rw6_n_h2 = feols(lhs_new ~ 1|year + race2 + my_geo + age +race2^my_geo + age^my_geo +year^my_geo + 
                        age^race2 + year^race2 +race2^my_geo^year + race2^my_geo^age  |
                        logrw ~ laglogrw, data = data_loc2, vcov = "hetero")
summary(fe_iv_rw6_n_h2)

fe_iv_rw7_n2 = feols(lhs_new ~ 1|year + race2 + my_geo + age +race2^my_geo + age^my_geo +year^my_geo + 
                      age^race2 + year^race2  + race2^my_geo^age  |
                      logrw ~ laglogrw, data = data_loc2, cluster = c("my_geo"))
summary(fe_iv_rw7_n2)

fe_iv_rw7_n_h2 = feols(lhs_new ~ 1|year + race2 + my_geo + age +race2^my_geo + age^my_geo +year^my_geo + 
                         age^race2 + year^race2  + race2^my_geo^age  |
                         logrw ~ laglogrw, data = data_loc2, vcov = "hetero")
summary(fe_iv_rw7_n_h2)

etable(fe_iv_rw3_n2, fe_iv_rw7_n2, fe_iv_rw6_n2)

# suppose that the migration elasticity 1/nu is 1.310
nu = 1/1.319


# using nu induced by the estimation, get the migration costs
mig_cost_df2 = mutate(mig_cost_df2, tau = -nu * tau_tilde)

# Julia requires UTF-8 encoding
write.csv(mig_cost_df2, "output/csv/migration_cost_oridess.csv", fileEncoding = "UTF-8")

# back out amenities
# primarily use the result of fe_iv_rw6_n2
fe_iv_rw6_n2_fe = fixef(fe_iv_rw6_n2)
fe_iv_rw6_n2_fe$year
fe_iv_rw6_n2_fe$race2
fe_iv_rw6_n2_fe$my_geo
fe_iv_rw6_n2_fe$age
fe_iv_rw6_n2_fe$`race2^my_geo`
fe_iv_rw6_n2_fe$`age^my_geo`
fe_iv_rw6_n2_fe$`year^my_geo`
fe_iv_rw6_n2_fe$`age^race2`
fe_iv_rw6_n2_fe$`year^race2`
fe_iv_rw6_n2_fe$`race2^my_geo^year`
fe_iv_rw6_n2_fe$`race2^my_geo^age`

# dataframe to store amenities (or estimated fixed effects in the second step)
ame_df = select(data2, race2, year, age, my_geo)
ame_df = unique(ame_df)
ame_df = ame_df[order(ame_df$race2, ame_df$year, ame_df$age, ame_df$my_geo), ]

# year
year_fe = data.frame(as.integer(names(fe_iv_rw6_n2_fe$year)), fe_iv_rw6_n2_fe$year)
names(year_fe) = c("year", "year_fe")

# race2
race2_fe = data.frame(as.integer(names(fe_iv_rw6_n2_fe$race2)), fe_iv_rw6_n2_fe$race2)
names(race2_fe) = c("race2", "race2_fe")

# my_geo
my_geo_fe = data.frame(as.integer(names(fe_iv_rw6_n2_fe$my_geo)), fe_iv_rw6_n2_fe$my_geo)
names(my_geo_fe) = c("my_geo", "my_geo_fe")

# age
age_fe = data.frame(as.integer(names(fe_iv_rw6_n2_fe$age)), fe_iv_rw6_n2_fe$age)
names(age_fe) = c("age", "age_fe")

# race2^my_geo
race2_my_geo_fe = data.frame(names(fe_iv_rw6_n2_fe$`race2^my_geo`), fe_iv_rw6_n2_fe$`race2^my_geo`)
names(race2_my_geo_fe) = c("race2_my_geo", "race2_my_geo_fe")
split_race2_my_geo = str_split(race2_my_geo_fe$race2_my_geo , "_", simplify = TRUE)
split_race2_my_geo = as.data.frame(split_race2_my_geo)

names(split_race2_my_geo) = c("race2", "my_geo")
split_race2_my_geo$race2 = as.integer(split_race2_my_geo$race2)
split_race2_my_geo$my_geo = as.integer(split_race2_my_geo$my_geo)

race2_my_geo_fe = cbind(race2_my_geo_fe, split_race2_my_geo)
race2_my_geo_fe = select(race2_my_geo_fe, -race2_my_geo)

rm(split_race2_my_geo)

# age^my_geo
age_my_geo_fe = data.frame(names(fe_iv_rw6_n2_fe$`age^my_geo`), fe_iv_rw6_n2_fe$`age^my_geo`)
names(age_my_geo_fe) = c("age_my_geo", "age_my_geo_fe")
split_age_my_geo = str_split(age_my_geo_fe$age_my_geo , "_", simplify = TRUE)
split_age_my_geo = as.data.frame(split_age_my_geo)

names(split_age_my_geo) = c("age", "my_geo")
split_age_my_geo$age = as.integer(split_age_my_geo$age)
split_age_my_geo$my_geo = as.integer(split_age_my_geo$my_geo)

age_my_geo_fe = cbind(age_my_geo_fe, split_age_my_geo)
age_my_geo_fe = select(age_my_geo_fe, -age_my_geo)

rm(split_age_my_geo)

# year^my_geo
year_my_geo_fe = data.frame(names(fe_iv_rw6_n2_fe$`year^my_geo`), fe_iv_rw6_n2_fe$`year^my_geo`)
names(year_my_geo_fe) = c("year_my_geo", "year_my_geo_fe")
split_year_my_geo = str_split(year_my_geo_fe$year_my_geo , "_", simplify = TRUE)
split_year_my_geo = as.data.frame(split_year_my_geo)

names(split_year_my_geo) = c("year", "my_geo")
split_year_my_geo$year = as.integer(split_year_my_geo$year)
split_year_my_geo$my_geo = as.integer(split_year_my_geo$my_geo)

year_my_geo_fe = cbind(year_my_geo_fe, split_year_my_geo)
year_my_geo_fe = select(year_my_geo_fe, -year_my_geo)

rm(split_year_my_geo)

# age^race2
age_race2_fe = data.frame(names(fe_iv_rw6_n2_fe$`age^race2`), fe_iv_rw6_n2_fe$`age^race2`)
names(age_race2_fe) = c("age_race2", "age_race2_fe")
split_age_race2 = str_split(age_race2_fe$age_race2 , "_", simplify = TRUE)
split_age_race2 = as.data.frame(split_age_race2)

names(split_age_race2) = c("age", "race2")
split_age_race2$age = as.integer(split_age_race2$age)
split_age_race2$race2 = as.integer(split_age_race2$race2)

age_race2_fe = cbind(age_race2_fe, split_age_race2)
age_race2_fe = select(age_race2_fe, -age_race2)

rm(split_age_race2)

# year^race2
year_race2_fe = data.frame(names(fe_iv_rw6_n2_fe$`year^race2`), fe_iv_rw6_n2_fe$`year^race2`)
names(year_race2_fe) = c("year_race2", "year_race2_fe")
split_year_race2 = str_split(year_race2_fe$year_race2 , "_", simplify = TRUE)
split_year_race2 = as.data.frame(split_year_race2)

names(split_year_race2) = c("year", "race2")
split_year_race2$year = as.integer(split_year_race2$year)
split_year_race2$race2 = as.integer(split_year_race2$race2)

year_race2_fe = cbind(year_race2_fe, split_year_race2)
year_race2_fe = select(year_race2_fe, -year_race2)

rm(split_year_race2)

# race2^my_geo^year
race2_my_geo_year_fe = data.frame(names(fe_iv_rw6_n2_fe$`race2^my_geo^year`), fe_iv_rw6_n2_fe$`race2^my_geo^year`)
names(race2_my_geo_year_fe) = c("race2_my_geo_year", "race2_my_geo_year_fe")
split_race2_my_geo_year = str_split(race2_my_geo_year_fe$race2_my_geo_year , "_", simplify = TRUE)
split_race2_my_geo_year = as.data.frame(split_race2_my_geo_year)

names(split_race2_my_geo_year) = c("race2", "my_geo", "year")
split_race2_my_geo_year$race2 = as.integer(split_race2_my_geo_year$race2)
split_race2_my_geo_year$my_geo = as.integer(split_race2_my_geo_year$my_geo)
split_race2_my_geo_year$year = as.integer(split_race2_my_geo_year$year)

race2_my_geo_year_fe = cbind(race2_my_geo_year_fe, split_race2_my_geo_year)
race2_my_geo_year_fe = select(race2_my_geo_year_fe, -race2_my_geo_year)

rm(split_race2_my_geo_year)

# race2^my_geo^age
race2_my_geo_age_fe = data.frame(names(fe_iv_rw6_n2_fe$`race2^my_geo^age`), fe_iv_rw6_n2_fe$`race2^my_geo^age`)
names(race2_my_geo_age_fe) = c("race2_my_geo_age", "race2_my_geo_age_fe")
split_race2_my_geo_age = str_split(race2_my_geo_age_fe$race2_my_geo_age , "_", simplify = TRUE)
split_race2_my_geo_age = as.data.frame(split_race2_my_geo_age)

names(split_race2_my_geo_age) = c("race2", "my_geo", "age")
split_race2_my_geo_age$race2 = as.integer(split_race2_my_geo_age$race2)
split_race2_my_geo_age$my_geo = as.integer(split_race2_my_geo_age$my_geo)
split_race2_my_geo_age$age = as.integer(split_race2_my_geo_age$age)

race2_my_geo_age_fe = cbind(race2_my_geo_age_fe, split_race2_my_geo_age)
race2_my_geo_age_fe = select(race2_my_geo_age_fe, -race2_my_geo_age)

rm(split_race2_my_geo_age)

# collect fixed effects in a single dataframe
ame_df = left_join(ame_df, year_fe, by = "year")
ame_df = left_join(ame_df, race2_fe, by = "race2")
ame_df = left_join(ame_df, my_geo_fe, by = "my_geo")
ame_df = left_join(ame_df, age_fe, by = "age")
ame_df = left_join(ame_df, race2_my_geo_fe, by = c("race2", "my_geo"))
ame_df = left_join(ame_df, age_my_geo_fe, by = c("age", "my_geo"))
ame_df = left_join(ame_df, year_my_geo_fe, by = c("year", "my_geo"))
ame_df = left_join(ame_df, age_race2_fe, by = c("age", "race2"))
ame_df = left_join(ame_df, year_race2_fe, by = c("year", "race2"))
ame_df = left_join(ame_df, race2_my_geo_year_fe, by = c("race2", "my_geo", "year"))
ame_df = left_join(ame_df, race2_my_geo_age_fe, by = c("race2", "my_geo", "age"))

ame_df = mutate(ame_df, B_tilde = year_fe + race2_fe + my_geo_fe + age_fe +
                  race2_my_geo_fe + age_my_geo_fe + year_my_geo_fe + age_race2_fe +
                  year_race2_fe + race2_my_geo_year_fe + race2_my_geo_age_fe)

ame_df = mutate(ame_df,
                B = exp(nu*B_tilde))

# amenities that are not normalized
# ame_nm_df = ame_df
# 
# write.csv(ame_nm_df, "output/csv/amenities_nm.csv", fileEncoding = "UTF-8")

ame_df = ame_df %>%
  group_by(race2, age, year) %>%
  mutate(meanB = mean(B)) %>%
  mutate(normalizedB = B / meanB)

write.csv(ame_df, "output/csv/amenities_oridess.csv", fileEncoding = "UTF-8")

# I can compute amenities for ages 20-70, years 1950-2000 in this way
ame_df_val = filter(ame_df, age %in% c(20, 30, 40, 50, 60, 70),
                    year %in% c(1950, 1960, 1970, 1980, 1990, 2000))

write.csv(ame_df_val, "output/csv/amenities_val_oridess.csv", fileEncoding = "UTF-8")


# ame_nm_df_val = filter(ame_nm_df, age %in% c(20, 30, 40, 50, 60, 70),
#                     year %in% c(1950, 1960, 1970, 1980, 1990, 2000))
# 
# write.csv(ame_nm_df_val, "output/csv/amenities_nm_val.csv", fileEncoding = "UTF-8")
