## Read and harmonise the raw register sources (Border Guard, Police, prison, ZUS, PESEL)
## into one linked table. Produces `full_database` and `full_database_processed`
## (18+ non-Schengen modelling frame) and writes the cleaned data/*.csv files.

## codes for schengen
schengen_iso3 <- c( "AUT", "BEL", "BGR", "HRV", "CZE", "DNK", "EST", "FIN", "FRA", "DEU", "GRC", "HUN", "ISL", "ITA", "LVA", "LIE", "LTU", "LUX", "MLT", "NLD", "NOR", "POL", "PRT", "ROU", "SVK", "SVN", "ESP", "SWE", "CHE")


## Polish border guards data -- within country
border <- fread("data-raw/undocumented_stay_2014_2024_country.csv")
setnames(border, c("citizenship", "count"), c("country_code", "border"))

border[, country:=countrycode(country_code, "iso3c", "country.name")]
border[country_code %in% c("RKS", "KOSOVO"), country:="KOSOVO"]
border[country_code == "XXX", country:="Stateless"]
border[country_code == "UNK", country:="unknown"]
border[, continent := countrycode(country_code, "iso3c", "continent")]
border[country == "KOSOVO", continent := "Europe"]
border[country_code == "IOT", continent := "Africa"]  #  British Indian Ocean Territory
border[is.na(continent), continent := "other/unknown"]
border[, wbd_region7 := countrycode(country_code, "iso3c", "region")]
border[country == "KOSOVO", wbd_region7 := "Europe & Central Asia"]
border[country_code == "IOT", wbd_region7 := "Sub-Saharan Africa"]   # matches MUS/SYC
border[is.na(wbd_region7), wbd_region7 := "other/unknown"]
border[, schengen := ifelse(country_code %in% schengen_iso3, "Schengen", "non-Schengen")]
border[country_code == "HRV" & year <  2023 & schengen == "Schengen", schengen := "non-Schengen"]
border[country_code %in% c("BGR","ROU") & year <  2025 & schengen == "Schengen", schengen := "non-Schengen"]

border[, sex := stri_trans_toupper(stri_split_fixed(sex_age, "_", n=2, simplify=T)[, 1])]
border[sex == "K", sex := "F"]
border[, sex:=tolower(sex)]
border[, age := stri_split_fixed(sex_age, "_", n=2, simplify=T)[, 2]]
border[, sex_age := NULL]
border[country_code == "PRK", ":="(country_code="KOR", country="South Korea")]


# Polish border guards data -- Information on borders

border_country <- fread("data-raw/undocumented_stay_2014_2024_borders.csv")
setnames(border_country, c("citizenship", "count"), c("country_code", "border"))

border_country[, country:=countrycode(country_code, "iso3c", "country.name")]
border_country[country_code %in% c("RKS", "KOSOVO"), country:="KOSOVO"]
border_country[country_code == "XXX", country:="Stateless"]
border_country[country_code == "UNK", country:="unknown"]
border_country[, continent := countrycode(country_code, "iso3c", "continent")]
border_country[country == "KOSOVO", continent := "Europe"]
border_country[country_code == "IOT", continent := "Africa"]  #  British Indian Ocean Territory
border_country[is.na(continent), continent := "other/unknown"]
border_country[, wbd_region7 := countrycode(country_code, "iso3c", "region")]
border_country[country == "KOSOVO", wbd_region7 := "Europe & Central Asia"]
border_country[country_code == "IOT", wbd_region7 := "Sub-Saharan Africa"]   # matches MUS/SYC
border_country[is.na(wbd_region7), wbd_region7 := "other/unknown"]
border_country[, schengen := ifelse(country_code %in% schengen_iso3, "Schengen", "non-Schengen")]
border_country[country_code == "HRV" & year <  2023 & schengen == "Schengen", schengen := "non-Schengen"]
border_country[country_code %in% c("BGR","ROU") & year <  2025 & schengen == "Schengen", schengen := "non-Schengen"]

border_country[, sex := stri_trans_toupper(stri_split_fixed(sex_age, "_", n=2, simplify=T)[, 1])]
border_country[sex == "K", sex := "F"]
border_country[, sex:=tolower(sex)]
border_country[, age := stri_split_fixed(sex_age, "_", n=2, simplify=T)[, 2]]
border_country[, sex_age := NULL]
border_country[country_code == "PRK", ":="(country_code="KOR", country="South Korea")]

## Police data

## NOTE (public replication): the raw police extract is individual, proceeding-level
## micro-data (one row per person/incident, with a unique proceeding id, and including
## criminal-proceeding records). It is a restricted source and is NOT shared here. Only the
## aggregated public table `data/poland-police.csv` (counts by year/country/age/sex, no id)
## ships, so we load that pre-built aggregate instead of rebuilding it from the raw file.
## The raw file `data-raw/police-data-2016-2024.csv` and the block that aggregated it
## (subset to criminal/process.reg, dedup by id, count by PESEL flag) are omitted. The
## aggregation logic is preserved in the full (non-public) package. See README.md.
police_data <- fread("data/poland-police.csv")

## prison data

#> prison[,.N, category]
#    category     N
#      <char> <int>
#1: temporary   616
#2: convicted   573
#3:  punished    52

prison <- fread("data-raw/prison-data-2010-2025.csv") |>
  subset(year %in% 2016:2025 & month == 12)

prison_data <- prison[, .(prison = sum(count)), .(year, sex, citizenship)]
prison_data[nchar(citizenship) ==0, citizenship := "XXX"]
setnames(prison_data, "citizenship", "country_code")

prison_data[, country:=countrycode(country_code, "iso3c", "country.name")]
prison_data[country_code == "XXX", country:="Stateless"]
prison_data[country_code == "XKX", country:="KOSOVO"]
prison_data[, continent := countrycode(country_code, "iso3c", "continent")]
prison_data[country == "KOSOVO", continent := "Europe"]
prison_data[is.na(continent), continent := "other/unknown"]
prison_data[, wbd_region7 := countrycode(country_code, "iso3c", "region")]
prison_data[country == "KOSOVO", wbd_region7 := "Europe & Central Asia"]
prison_data[is.na(wbd_region7), wbd_region7 := "other/unknown"]
prison_data[, schengen := ifelse(country_code %in% schengen_iso3, "Schengen", "non-Schengen")]
prison_data[country_code == "HRV" & year <  2023 & schengen == "Schengen", schengen := "non-Schengen"]
prison_data[country_code %in% c("BGR","ROU") & year <  2025 & schengen == "Schengen", schengen := "non-Schengen"]


## social insurance data
zus <- fread("data-raw/social-insurance_2016_2022q3_2024.csv")
## 2016-2021 AND 2024 share the same definition (disability.pension / all.people, Q4) in this file;
## zus_2024.csv was a different measure (health insurance) so the big file is used for 2024 instead.
zus_16_21_24 <- subset(zus,
                    type == "disability.pension" & group == "all.people" &
                    quarter == 4 & year %in% c(2016:2021, 2024) &
                    !age %in% c("17", "19"))

zus2022 <- fread("data-raw/social-insurance_2022.csv")
zus2023 <- fread("data-raw/social-insurance_2023.csv")
zus_22_23 <- rbind(zus2022, zus2023)[age >= 18 & age <= 100]   # cap implausible ages

zus_data <- rbind(
  zus_16_21_24[, .(zus = sum(count, na.rm=T)), keyby=.(year, country_code, sex)],
  zus_22_23[, .(zus = sum(count)), keyby=.(year, country_code=citizenship, sex)]
)

zus_data[country_code == "SCG", country_code := "SRB"]   # collapse Serbia & Montenegro into Serbia

zus_data[, country:=countrycode(country_code, "iso3c", "country.name")]
zus_data[country_code == "XXX", country:="unknown"]       # ZUS 'other' category (2016-2021, 2024)
zus_data[country_code == "UNKNOWN", country:="unknown"]   # same category, 2022-2023 files
zus_data[country_code == "YUG", country:="Yugoslavia"]

zus_data[, continent := countrycode(country_code, "iso3c", "continent")]
zus_data[country_code %in% c("IOT","ATF"), continent := "Africa"]  # British Indian Ocean & French Southern Territories
zus_data[country %in% c("KOSOVO", "Yugoslavia", "Soviet Union", "Czechoslovakia"), continent := "Europe"]
zus_data[is.na(continent), continent := "other/unknown"]
zus_data[, wbd_region7 := countrycode(country_code, "iso3c", "region")]
zus_data[country_code %in% c("IOT","ATF"), wbd_region7 := "Sub-Saharan Africa"]   # Indian Ocean; matches MUS/SYC/REU
zus_data[country %in% c("KOSOVO", "Yugoslavia", "Soviet Union", "Czechoslovakia",
                        "Svalbard & Jan Mayen"),
         wbd_region7 := "Europe & Central Asia"]
zus_data[country %in% c("Wallis & Futuna"), wbd_region7 := "East Asia & Pacific"]
zus_data[is.na(wbd_region7), wbd_region7 := "other/unknown"]
zus_data[, schengen := ifelse(country_code %in% schengen_iso3, "Schengen", "non-Schengen")]
zus_data[country_code == "HRV" & year <  2023 & schengen == "Schengen", schengen := "non-Schengen"]
zus_data[country_code %in% c("BGR","ROU") & year <  2025 & schengen == "Schengen", schengen := "non-Schengen"]

## koreas are together
zus_data[country_code == "PRK", ":="(country_code="KOR", country="South Korea")]

# PESEL register
pesel_17_22 <- fread("data-raw/pesel-register-2017-2022.csv")
pesel_17_22 <- subset(pesel_17_22,
       year %in% 2017:2021 & quarter == 4 &
       !age %in% c("to 17", "17"))
pesel_22 <- fread("data-raw/pesel-register-2022-2023.csv")[age >= 18 & age <= 100]
pesel_24 <- fread("data-raw/pesel-register-2024.csv")[country != ""]

pesel_24[pesel_17_22[, .N, keyby=.(country, iso3c_grok)], on = "country", country_code := i.iso3c_grok]
pesel_24[country == "MIKRONEZYJSKIE", country_code:="FSM"]
pesel_24[country == "SALOMOŃSKIE", country_code:="SLB"]
pesel_24[country == "WATYKAŃSKIE", country_code:="VAT"]
pesel_24[country == "WYSPY MINOR", country_code:="UMI"]
pesel_24[country == "ANTIGUAŃSKO-BARBUDZKIE", country_code:="ATG"]
## historical-citizenship names not present in the 2017-2022 lookup
pesel_24[country == "CZECHOSŁOWACKIE", country_code:="CSK"]
pesel_24[country == "JUGOSŁOWIAŃSKIE", country_code:="YUG"]
pesel_24[country == "NIEMIECKIE/BERLIN ZACH.", country_code:="DEU"]
pesel_24[country == "NIEMIECKIE/NRD", country_code:="DDR"]
pesel_24[country == "RADZIECKIE/ZSRR", country_code:="SUN"]
pesel_24[is.na(country_code), country_code := "UNKNOWN"]

pesel_data <- rbind(
  pesel_17_22[, .(pesel=sum(pesel)), .(year, country, country_code=iso3c_grok, sex=ifelse(sex %in% c("1", "m"), "m", "f"))],
  pesel_22[, .(year, country_code=citizenship, sex, pesel=count)],
  pesel_24[, .(year=2024, country_code, sex, pesel)],
  fill = T
)

pesel_data[nchar(country_code) == 0, country_code := "ATG"]   # only blank-code name is Antigua & Barbuda

pesel_data[, country:=countrycode(country_code, "iso3c", "country.name")]
pesel_data[country_code %in% c("XKX","KOSOVO"), country:="KOSOVO"]
pesel_data[country_code == "XXA", country:="Stateless"]       # PESEL: XXA = stateless
pesel_data[country_code %in% c("XXX","UNKNOWN"), country:="unknown"]  # PESEL: XXX = other/unknown
pesel_data[country_code == "YUG", country:="Yugoslavia"]
pesel_data[country_code == "SCG", country:="Serbia and Montenegro"]
pesel_data[country_code == "SUN", country:="Soviet Union"]
pesel_data[country_code == "CSK", country:="Czechoslovakia"]
pesel_data[country_code == "DDR", country:="East Germany"]

pesel_data[, continent := countrycode(country_code, "iso3c", "continent")]
pesel_data[country_code %in% c("IOT","ATF"), continent := "Africa"]  # British Indian Ocean & French Southern Territories
pesel_data[country_code == "UMI", continent := "Oceania"]  # US Minor Outlying Islands
pesel_data[country %in% c("KOSOVO","Yugoslavia","Serbia and Montenegro","Soviet Union","Czechoslovakia","East Germany"), continent := "Europe"]
pesel_data[is.na(continent), continent := "other/unknown"]
pesel_data[, wbd_region7 := countrycode(country_code, "iso3c", "region")]
pesel_data[country_code %in% c("IOT","ATF"), wbd_region7 := "Sub-Saharan Africa"]   # Indian Ocean; matches MUS/SYC/REU
pesel_data[country_code == "UMI", wbd_region7 := "East Asia & Pacific"]
pesel_data[country %in% c("KOSOVO","Yugoslavia","Serbia and Montenegro","Soviet Union","Czechoslovakia","East Germany"), wbd_region7 := "Europe & Central Asia"]
pesel_data[is.na(wbd_region7), wbd_region7 := "other/unknown"]
pesel_data[, schengen := ifelse(country_code %in% schengen_iso3, "Schengen", "non-Schengen")]
pesel_data[country_code == "HRV" & year <  2023 & schengen == "Schengen", schengen := "non-Schengen"]
pesel_data[country_code %in% c("BGR","ROU") & year <  2025 & schengen == "Schengen", schengen := "non-Schengen"]
pesel_data[country_code == "PRK", ":="(country_code="KOR", country="South Korea")]


## NOTE (public replication): the tax register is a restricted source and is NOT shared here.
## The whole tax block is therefore omitted; every tax-based output is dropped downstream
## (Figure 2 Tax panel, Table 2 Tax column, appendix Tax register diagnostics, and the Tax
## rows of the cross-register AIC/popsize comparisons). See README.md in this folder.



## Linking all datasets into one

codes_unique <- unique(c(border$country_code,
                         police_data$country_code,
                         prison_data$country_code,
                         zus_data$country_code,
                         pesel_data$country_code))

years_unique <- unique(c(border$year,
                         police_data$year,
                         prison_data$year,
                         zus_data$year,
                         pesel_data$year))

sex_unique <- unique(c(border$sex,
                         police_data$sex,
                         prison_data$sex,
                         zus_data$sex,
                         pesel_data$sex))

full_database <- CJ(year=years_unique, sex=sex_unique, country_code=codes_unique)

full_database[border[same_year == FALSE & age != "17", .(border=sum(border)), .(year, country_code, sex)],
              on = c("year", "sex", "country_code"),
              border := i.border]

full_database[police_data[, .(police_id_no=sum(police_id_no), 
                              police_id_yes=sum(police_id_yes),
                              police = sum(police_id_yes+police_id_no)), .(year, country_code, sex)],
              on = c("year", "sex", "country_code"),
              ":="(police_id_no=i.police_id_no, police_id_yes=i.police_id_yes, police=i.police)]

full_database[prison_data[, .(prison=sum(prison)), .(year, sex, country_code)],
              on = c("year", "sex", "country_code"),
              prison := i.prison]

full_database[zus_data[, .(zus=sum(zus)), .(year, sex, country_code)],
              on = c("year", "sex", "country_code"),
              pop_insured := i.zus]

full_database[pesel_data[, .(pesel=sum(pesel)), .(year, sex, country_code)],
              on = c("year", "sex", "country_code"),
              pop_register := i.pesel]

## pop_tax join omitted -- tax data restricted (public replication)

full_database[, lapply(.SD, \(x) sum(is.na(x))), keyby=.(year,sex), .SDcols = border:pop_register]

## harmonise country codes across sources
full_database[country_code == "SCG", country_code := "SRB"]
full_database[, country:=countrycode(country_code, "iso3c", "country.name")]
full_database[, country_code2 := country_code]
full_database[country_code %in% c("XXX", "SOT", "ANT", "UNKNOWN", "UNK", "XXA", "DDR", "YUG", "CSK", "SUN"), country_code2 := "unknown/stateless"]
full_database[country_code %in% c("XXX", "SOT", "ANT", "UNKNOWN", "UNK", "XXA", "DDR", "YUG", "CSK", "SUN"), country := "unknown/stateless"]
full_database[country_code %in% c("RKS", "KOSOVO", "XKX"), country:="KOSOVO"]
full_database[country_code %in% c("RKS", "KOSOVO", "XKX"), country_code2:="KOSOVO"]
full_database[country_code == "PRK", ":="(country_code2="KOR", country="South Korea")]
full_database[, schengen := ifelse(country_code %in% schengen_iso3, "Schengen", "non-Schengen")]
full_database[country_code == "HRV" & year <  2023 & schengen == "Schengen", schengen := "non-Schengen"]
full_database[country_code %in% c("BGR","ROU") & year <  2025 & schengen == "Schengen", schengen := "non-Schengen"]


## aggregate to harmonised codes and classify for the model
full_database_processed <- full_database[year %in% 2019:2024, lapply(.SD, sum, na.rm=T),
                                         .(year, sex, country_code=country_code2, country), .SDcols = border:pop_register]

## Drop EU free-movement nationals from the modelling/table population: here "non-Schengen"
## means non-EU (genuine third countries). EU citizens have free movement and are absent
## from the apprehension data (m = 0) while carrying large reference populations. The rule
## covers HRV/BGR/ROU and the EU members outside Schengen, IRL and CYP, throughout, and
## GBR up to the end of the Brexit transition period (free movement through 2020; a
## visa-free third country from 2021). The schengen classification is kept intact in
## `full_database`; this exclusion (noted in the paper) only applies to
## `full_database_processed` and everything built from it.
EU_FREE_MOVEMENT <- c("HRV", "BGR", "ROU", "IRL", "CYP")
full_database_processed <- full_database_processed[!country_code %in% EU_FREE_MOVEMENT]
full_database_processed <- full_database_processed[!(country_code == "GBR" & year <= 2020)]

full_database_processed[, continent := countrycode(country_code, "iso3c", "continent")]
full_database_processed[, wbd_region7 := countrycode(country_code, "iso3c", "region")]
full_database_processed[, schengen := ifelse(country_code %in% schengen_iso3, "Schengen", "non-Schengen")]
full_database_processed[country_code == "HRV" & year <  2023 & schengen == "Schengen", schengen := "non-Schengen"]
full_database_processed[country_code %in% c("BGR","ROU") & year <  2025 & schengen == "Schengen", schengen := "non-Schengen"]

full_database_processed[country_code == "KOSOVO", ":="(continent="Europe",  wbd_region7="Europe & Central Asia")]
full_database_processed[country_code %in% c("ATF","IOT"), ":="(continent="Africa",  wbd_region7="Sub-Saharan Africa")]
full_database_processed[country_code == "UMI", ":="(continent="Oceania", wbd_region7="East Asia & Pacific")]
full_database_processed[country_code == "ESH", wbd_region7 := "Middle East & North Africa"]   # Western Sahara
full_database_processed[country_code == "PRI", wbd_region7 := "Latin America & Caribbean"]      # Puerto Rico
full_database_processed[is.na(continent), ":="(continent ="unknown/stateless", wbd_region7 ="unknown/stateless")]
full_database_processed[is.na(wbd_region7), wbd_region7 := "unknown/stateless"]   # e.g. Antarctica: no WB region
full_database_processed[, ukr := fifelse(country_code == "UKR", 1L, 0L)]
full_database_processed[, simplified_proc := fifelse(
  country_code %in% c("UKR", "BLR", "MDA", "GEO", "ARM", "RUS"), 1L, 0L)]

full_database_processed <- full_database_processed[schengen == "non-Schengen"]


## write all
fwrite(border, file = "data/poland-unauthorized.csv")
fwrite(border_country, file = "data/poland-unauthorized-borders.csv")
## poland-police.csv is shipped as a pre-built public aggregate (the raw police micro-data is
## restricted and omitted), so it is read above as an input and not re-written here.
fwrite(prison_data, file = "data/poland-prisons.csv")
fwrite(pesel_data, file = "data/poland-population-register.csv")
fwrite(zus_data, file = "data/poland-social-insurance.csv")
fwrite(full_database, file = "data/poland-full-database.csv")
fwrite(full_database_processed, file = "data/poland-for-model.csv")

