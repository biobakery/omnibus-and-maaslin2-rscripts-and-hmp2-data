# Omnibus tests on mgx and mtx data from the HMP2. Baseline only data.

dir.create("R_Metagenomic_Statistics") # Create a new directory
dir.create("R_Metagenomic_Statistics/Data") # Create a new directory
setwd("R_Metagenomic_Statistics") # Change the current working directory 

# Load the packages needed
library(vegan)
library(plyr)

#####
# MGX species
#####
# You can either download the file from bitbucket page or do it with the below code:
download.file("https://github.com/biobakery/omnibus-and-maaslin2-rscripts-and-hmp2-data/blob/master/taxonomic_profiles_pcl_week0.csv", "./Data/taxonomic_profiles_pcl_week0.csv") # Download the mgx species data and put it into the data directory

# Read the taxonomic data into R environment
tax = read.csv(file = "./Data/taxonomic_profiles_pcl_week0.csv", header = T, row.names = 1, check.names = FALSE)

head(tax) # Inspect the tax
dim(tax) # dimensions of the tax
str(tax) # structure of the tax
names(tax) # column names of tax
row.names(tax) # row.names of tax


# Prepare the data
# Extract the metadata
metadata = data.frame(tax[1:5])
metadata[1:5,] # check the output
str(metadata) #structure of metadata
is.na(metadata) # Check for NAs that will mess with the PERMANOVAS: Age has some
count(is.na(metadata$consent_age)) # Check for how many there are: 6/96
# If this was a discrete variable we could just classify the NAs as Unknown and keep them in the model, 
# but since Age is a continuous variable typically we would either remove those from the data or impute the median. 
# In this case let's impute the median in order to keep samples. 
unique(metadata$consent_age)
metadata$consent_age[is.na(metadata$consent_age)] = median(metadata$consent_age, na.rm = T)
unique(metadata$consent_age) # Check the output: good to go

# Extract species data and transpose the df
species = data.frame(t(tax[6:ncol(tax)]))
str(species) # everything is numeric, good to go
row.names(species)
species[1:8,1:4] # check the output
# subset to species only
# which don't have "t__"
tmp.ind = grep("\\|t__", rownames(species), invert = T) # grep the rows that do not include strain stratifications
tmp.ind # check the output
tmp = species[tmp.ind,]  # Create a new dataframe with only those row numbers
tmp.ind = grep("\\|s__", rownames(tmp)) # grep the rows that only include down to species stratifications
tmp.ind # check the output
species = tmp[tmp.ind,] # Create a new dataframe with only those row numbers
rm(tmp,tmp.ind) # remove temp files to clear up space
row.names(species) # Check the output to make sure that we only have species level stratifications

# trim species names
rownames(species) = gsub(".*\\|", "", rownames(species))
row.names(species) # Check the output, looks great
colSums(species) # Check the sample sums to make sure they are in proportion format (0-1) and are all ~1

# filter for beta div (we will keep species as is for alpha diversity)
dim(species)
dim(species[apply(species, 1, function(x) sum(x > 0.0001) > 0.1 * ncol(species)), ]) 
species_filt = species[apply(species, 1, function(x) sum(x > 0.0001) > 0.1 * ncol(species)), ]
#Let's transpose it for easier use downstream
species_filt = data.frame(t(species_filt), check.names = F)
species = data.frame(t(species), check.names = F)


bray_species = vegdist(species_filt, method = "bray") # Calculate Bray-Curtis dissimilarity


#####
# Alpha diversity

# use unfiltered data for alpha div
shannon = data.frame(diversity(species, index = "shannon"))
simpson = data.frame(diversity(species, index = "simpson"))
invsimpson = data.frame(diversity(species, index = "invsimpson"))
# merge all into one dataframe
alphadiv = cbind(shannon, simpson, invsimpson)
head(alphadiv) # check the output
# rename the columns to shorter names
names(alphadiv) = c("Shannon", "Simpson", "InvSimpson")
head(alphadiv) # check the output
# append metadata
alpha_met_df <- merge(alphadiv, metadata, by = "row.names")
row.names(alpha_met_df) = alpha_met_df$Row.names # Make the samples the rownames again
alpha_met_df$Row.names = NULL # Get rid of the Row.names column
head(alpha_met_df)

### Shannon alpha diversity
#### Univariate

# Can do this with a for loop
for (col in names(alpha_met_df)[-c(2:3)]){
  alpha_shannon_univ = anova(lm(as.formula(paste("Shannon ~ ", col)), data = alpha_met_df[-c(2:3)])) 
  print(alpha_shannon_univ)
}

# Alternatively, can do it one-by-one
shannon_site = anova(lm(Shannon ~ site_name, data = alpha_met_df))
shannon_site

shannon_sex = anova(lm(Shannon ~ sex, data = alpha_met_df))
shannon_sex

shannon_race = anova(lm(Shannon ~ race, data = alpha_met_df))
shannon_race

shannon_age = anova(lm(Shannon ~ consent_age, data = alpha_met_df))
shannon_age

shannon_diagnosis = anova(lm(Shannon ~ diagnosis, data = alpha_met_df))
shannon_diagnosis

#### Multivariate
# Can do it without being verbose
shannon_multi = anova(lm(Shannon ~ ., data = alpha_met_df[-c(2:3)]))
shannon_multi

# Alternatively, can write out each variable in the model
shannon_multi = anova(lm(Shannon ~ site_name + sex + race + consent_age + diagnosis, data = alpha_met_df))
shannon_multi


### Simpson alpha diversity
#### Univariate

# Can do this with a for loop
for (col in names(alpha_met_df)[-c(1,3)]){
  alpha_simpson_univ = anova(lm(as.formula(paste("Simpson ~ ", col)), data = alpha_met_df[-c(1,3)])) 
  print(alpha_simpson_univ)
}

# Alternatively, can do it one-by-one
simpson_site = anova(lm(Simpson ~ site_name, data = alpha_met_df))
simpson_site

simpson_sex = anova(lm(Simpson ~ sex, data = alpha_met_df))
simpson_sex

simpson_race = anova(lm(Simpson ~ race, data = alpha_met_df))
simpson_race

simpson_age = anova(lm(Simpson ~ consent_age, data = alpha_met_df))
simpson_age

simpson_diagnosis = anova(lm(Simpson ~ diagnosis, data = alpha_met_df))
simpson_diagnosis

#### Multivariate
# Can do it without being verbose
simpson_multi = anova(lm(Simpson ~ ., data = alpha_met_df[-c(1,3)]))
simpson_multi

# Alternatively, can write out each variable in the model
simpson_multi = anova(lm(Simpson ~ site_name + sex + race + consent_age + diagnosis, data = alpha_met_df))
simpson_multi


### InvSimpson alpha diversity
#### Univariate

# Can do this with a for loop
for (col in names(alpha_met_df)[-c(1:2)]){
  alpha_invsimpson_univ = anova(lm(as.formula(paste("InvSimpson ~ ", col)), data = alpha_met_df[-c(1:2)])) 
  print(alpha_invsimpson_univ)
}

# Alternatively, can do it one-by-one
invsimpson_site = anova(lm(InvSimpson ~ site_name, data = alpha_met_df))
invsimpson_site

invsimpson_sex = anova(lm(InvSimpson ~ sex, data = alpha_met_df))
invsimpson_sex

invsimpson_race = anova(lm(InvSimpson ~ race, data = alpha_met_df))
invsimpson_race

invsimpson_age = anova(lm(InvSimpson ~ consent_age, data = alpha_met_df))
invsimpson_age

invsimpson_diagnosis = anova(lm(InvSimpson ~ diagnosis, data = alpha_met_df))
invsimpson_diagnosis


#### Multivariate
# Can do it without being verbose
invsimpson_multi = anova(lm(InvSimpson ~ ., data = alpha_met_df[-c(1:2)]))
invsimpson_multi

# Alternatively, can write out each variable in the model
invsimpson_multi = anova(lm(InvSimpson ~ site_name + sex + race + consent_age + diagnosis, data = alpha_met_df))
invsimpson_multi


# Now to look at all 3 measures together for diagnosis (univariate)
rbind(shannon_diagnosis, simpson_diagnosis, invsimpson_diagnosis)


##### 
#Beta Diversity: PERMANOVA Tests using Bray-Curtis Dissimilarity

### Univariate
# Can do this with a for loop
for (col in names(metadata)){
  adonis_univ <- adonis(as.formula(paste("bray_species ~ ", col)), data = metadata)
  print(adonis_univ)
}

# Alternatively, can do it one-by-one
adonis_site_tax = adonis(bray_species ~ site_name, data = metadata)
adonis_site_tax

adonis_sex_tax = adonis(bray_species ~ sex, data = metadata)
adonis_sex_tax

adonis_age_tax = adonis(bray_species ~ consent_age, data = metadata)
adonis_age_tax

adonis_race_tax = adonis(bray_species ~ race, data = metadata)
adonis_race_tax

adonis_diagnosis_tax = adonis(bray_species ~ diagnosis, data = metadata)
adonis_diagnosis_tax


#### Multivariate
# Can do it without being verbose
adonis_multi_tax = adonis(bray_species ~ ., data = metadata)
adonis_multi_tax

# Alternatively, can write out each variable in the model
adonis_multi_tax = adonis(bray_species ~ site_name + sex + race + consent_age + diagnosis, data = metadata)
adonis_multi_tax


#####
# MGX pathways
#####
# You can either download the file from bitbucket page or do it with the below code:
download.file("https://github.com/biobakery/omnibus-and-maaslin2-rscripts-and-hmp2-data/blob/master/dna_pathabundance_relab_pcl_week0.csv", "./Data/dna_pathabundance_relab_pcl_week0.csv")

# Read the dna pathway data into R environment
dna_path = read.csv(file = "./Data/dna_pathabundance_relab_pcl_week0.csv", header = T, row.names = 1, check.names = FALSE)

head(dna_path) # Inspect the data
dim(dna_path) # dimensions of the data
str(dna_path) # structure of the data
names(dna_path) # column names of data
row.names(dna_path) # row.names of data
# Remove metadata and keep only pathways and transpose the data
dna_path = data.frame(t(dna_path[6:ncol(dna_path)]))
str(dna_path) # everything is numeric, good to go
row.names(dna_path)

# Remove species stratifications
tmp.ind = grep("\\|.*", rownames(dna_path), invert = T) # grep the rows that do not include species stratifications 
tmp.ind # check the output
dna_path_unstratified = dna_path[tmp.ind,] # Create a new dataframe with only those unstratified rows
rm(tmp.ind) # Remove tmp.ind to clear space
row.names(dna_path_unstratified) # check the output: looks great
colSums(dna_path_unstratified) # Check the sample sums to make sure they are in proportion format (0-1) and are all ~1


# filter for beta div
dim(dna_path_unstratified)
dim(dna_path_unstratified[apply(dna_path_unstratified, 1, function(x) sum(x > 0.0001) > 0.1 * ncol(dna_path_unstratified)), ]) 
dna_path_unstratified_filt = dna_path_unstratified[apply(dna_path_unstratified, 1, function(x) sum(x > 0.0001) > 0.1 * ncol(dna_path_unstratified)), ]
#Let's transpose it for easier use downstream
dna_path_unstratified_filt = data.frame(t(dna_path_unstratified_filt), check.names = F)
dna_path_unstratified = data.frame(t(dna_path_unstratified), check.names = F)


bray_dna_path_unstratified = vegdist(dna_path_unstratified_filt, method = "bray") # Calculate Bray-Curtis dissimilarity

##### 
#Beta Diversity: PERMANOVA Tests using Bray-Curtis Dissimilarity

### Univariate
# Can do this with a for loop
for (col in names(metadata)){
  adonis_univ <- adonis(as.formula(paste("bray_dna_path_unstratified ~ ", col)), data = metadata)
  print(adonis_univ)
}

# Alternatively, can do it one-by-one
adonis_site_dna_path = adonis(bray_dna_path_unstratified ~ site_name, data = metadata)
adonis_site_dna_path

adonis_sex_dna_path = adonis(bray_dna_path_unstratified ~ sex, data = metadata)
adonis_sex_dna_path

adonis_age_dna_path = adonis(bray_dna_path_unstratified ~ consent_age, data = metadata)
adonis_age_dna_path

adonis_race_dna_path = adonis(bray_dna_path_unstratified ~ race, data = metadata)
adonis_race_dna_path

adonis_diagnosis_dna_path = adonis(bray_dna_path_unstratified ~ diagnosis, data = metadata)
adonis_diagnosis_dna_path


#### Multivariate
# Can do it without being verbose
adonis_multi_dna_path = adonis(bray_dna_path_unstratified ~ ., data = metadata)
adonis_multi_dna_path

# Alternatively, can write out each variable in the model
adonis_multi_dna_path = adonis(bray_dna_path_unstratified ~ site_name + sex + race + consent_age + diagnosis, data = metadata)
adonis_multi_dna_path


#####
# MTX pathways
#####
# You can either download the file from bitbucket page or do it with the below code:
download.file("https://github.com/biobakery/omnibus-and-maaslin2-rscripts-and-hmp2-data/blob/master/rna_pathabundance_relab_pcl_week0.csv", "./Data/rna_pathabundance_relab_pcl_week0.csv")

# Read the rna pathway data into R environment
rna_path = read.csv(file = "./Data/rna_pathabundance_relab_pcl_week0.csv", header = T, row.names = 1, check.names = FALSE)

head(rna_path) # Inspect the data
dim(rna_path) # dimensions of the data
str(rna_path) # structure of the data
names(rna_path) # column names of data
row.names(rna_path) # row.names of data
# Remove metadata and keep only pathways and transpose the data
rna_path = data.frame(t(rna_path[6:ncol(rna_path)]))
str(rna_path) # everything is numeric, good to go
row.names(rna_path)
# minimize the metadata to just the samples available in these data
dim(metadata)
list = names(rna_path) # make a list of sample ids to subset on
list # check the output
metadata_rna = subset(metadata, row.names(metadata) %in% list)
dim(metadata_rna)
metadata_rna # check the output

# Remove species stratifications
tmp.ind = grep("\\|.*", rownames(rna_path), invert = T) # grep the rows that do not include species stratifications 
tmp.ind # check the output
rna_path_unstratified = rna_path[tmp.ind,] # Create a new dataframe with only those unstratified rows
rm(tmp.ind) # Remove tmp.ind to clear space
row.names(rna_path_unstratified) # check the output: looks great
colSums(rna_path_unstratified) # Check the sample sums to make sure they are in proportion format (0-1) and are all ~1


# filter for beta div
dim(rna_path_unstratified)
dim(rna_path_unstratified[apply(rna_path_unstratified, 1, function(x) sum(x > 0.0001) > 0.1 * ncol(rna_path_unstratified)), ]) 
rna_path_unstratified_filt = rna_path_unstratified[apply(rna_path_unstratified, 1, function(x) sum(x > 0.0001) > 0.1 * ncol(rna_path_unstratified)), ]
#Let's transpose it for easier use downstream
rna_path_unstratified_filt = data.frame(t(rna_path_unstratified_filt), check.names = F)
rna_path_unstratified = data.frame(t(rna_path_unstratified), check.names = F)


bray_rna_path_unstratified = vegdist(rna_path_unstratified_filt, method = "bray") # Calculate Bray-Curtis dissimilarity

##### 
#Beta Diversity: PERMANOVA Tests using Bray-Curtis Dissimilarity

### Univariate
# Can do this with a for loop
for (col in names(metadata_rna)){
  adonis_univ <- adonis(as.formula(paste("bray_rna_path_unstratified ~ ", col)), data = metadata_rna)
  print(adonis_univ)
}

# Alternatively, can do it one-by-one
adonis_site_rna_path = adonis(bray_rna_path_unstratified ~ site_name, data = metadata_rna)
adonis_site_rna_path

adonis_sex_rna_path = adonis(bray_rna_path_unstratified ~ sex, data = metadata_rna)
adonis_sex_rna_path

adonis_age_rna_path = adonis(bray_rna_path_unstratified ~ consent_age, data = metadata_rna)
adonis_age_rna_path

adonis_race_rna_path = adonis(bray_rna_path_unstratified ~ race, data = metadata_rna)
adonis_race_rna_path

adonis_diagnosis_rna_path = adonis(bray_rna_path_unstratified ~ diagnosis, data = metadata_rna)
adonis_diagnosis_rna_path


#### Multivariate
# Can do it without being verbose
adonis_multi_rna_path = adonis(bray_rna_path_unstratified ~ ., data = metadata_rna)
adonis_multi_rna_path

# Alternatively, can write out each variable in the model
adonis_multi_rna_path = adonis(bray_rna_path_unstratified ~ site_name + sex + race + consent_age + diagnosis, data = metadata_rna)
adonis_multi_rna_path


#####
# RNA/DNA pathway ratios
#####
# You can either download the file from bitbucket page or do it with the below code:
download.file("https://github.com/biobakery/omnibus-and-maaslin2-rscripts-and-hmp2-data/blob/master/rna_dna_path_relative_expression_week0.csv", "./Data/rna_dna_path_relative_expression_week0.csv")

# Read the rna_dna pathway data into R environment
rna_dna_path = read.csv(file = "./Data/rna_dna_path_relative_expression_week0.csv", header = T, row.names = 1, check.names = FALSE)

head(rna_dna_path) # Inspect the data
dim(rna_dna_path) # dimensions of the data
str(rna_dna_path) # structure of the data
names(rna_dna_path) # column names of data
row.names(rna_dna_path) # row.names of data
# Transpose the data
rna_dna_path = data.frame(t(rna_dna_path))
str(rna_dna_path) # everything is numeric, good to go
row.names(rna_dna_path)
# minimize the metadata to just the samples available in these data
dim(metadata)
list = names(rna_dna_path) # make a list of sample ids to subset on
list # check the output
metadata_rna_dna = subset(metadata, row.names(metadata) %in% list)
dim(metadata_rna_dna)
metadata_rna_dna # check the output: we only have one race now, so let's get rid of that column
metadata_rna_dna$race = NULL

# Remove species stratifications
tmp.ind = grep("\\|.*", rownames(rna_dna_path), invert = T) # grep the rows that do not include species stratifications 
tmp.ind # check the output
rna_dna_path_unstratified = rna_dna_path[tmp.ind,] # Create a new dataframe with only those unstratified rows
rm(tmp.ind) # Remove tmp.ind to clear space
row.names(rna_dna_path_unstratified) # check the output: looks great


# filter for beta div
#Only keep RNA/DNA pathways that passed filtering for DNA pathways
#Create a list of the pathway names (col names) for subsetting the data
list = names(dna_path_unstratified_filt)
list
#Check the dimensions to make sure it matches with the DNA numbers before subsetting
dim(dna_path_unstratified_filt)
dim(subset(rna_dna_path_unstratified, row.names(rna_dna_path_unstratified) %in% list))
#subset
rna_dna_path_unstratified_filt = rna_dna_path_unstratified[list,]
dim(rna_dna_path_unstratified_filt)
head(rna_dna_path_unstratified_filt)

#Let's transpose the dataframe for easier use downstream
rna_dna_path_unstratified_filt = data.frame(t(rna_dna_path_unstratified_filt), check.names = F)

# log transform the RNA/DNA ratio
rna_dna_path_unstratified_filt_log = log2(rna_dna_path_unstratified_filt + 1)
head(rna_dna_path_unstratified_filt_log)


euclidean_rna_dna_path_unstratified = vegdist(rna_dna_path_unstratified_filt_log, method = "euclidean") # Calculate Euclidean distance on the pathway ratios

##### 
#Beta Diversity: PERMANOVA Tests using Euclidean distances

### Univariate
# Can do this with a for loop
for (col in names(metadata_rna_dna)){
  adonis_univ <- adonis(as.formula(paste("euclidean_rna_dna_path_unstratified ~ ", col)), data = metadata_rna_dna)
  print(adonis_univ)
}

# Alternatively, can do it one-by-one
adonis_site_rna_dna_path = adonis(euclidean_rna_dna_path_unstratified ~ site_name, data = metadata_rna_dna)
adonis_site_rna_dna_path

adonis_sex_rna_dna_path = adonis(euclidean_rna_dna_path_unstratified ~ sex, data = metadata_rna_dna)
adonis_sex_rna_dna_path

adonis_age_rna_dna_path = adonis(euclidean_rna_dna_path_unstratified ~ consent_age, data = metadata_rna_dna)
adonis_age_rna_dna_path

adonis_diagnosis_rna_dna_path = adonis(euclidean_rna_dna_path_unstratified ~ diagnosis, data = metadata_rna_dna)
adonis_diagnosis_rna_dna_path


#### Multivariate
# Can do it without being verbose
adonis_multi_rna_dna_path = adonis(euclidean_rna_dna_path_unstratified ~ ., data = metadata_rna_dna)
adonis_multi_rna_dna_path

# Alternatively, can write out each variable in the model
adonis_multi_rna_dna_path = adonis(euclidean_rna_dna_path_unstratified ~ site_name + sex + consent_age + diagnosis, data = metadata_rna_dna)
adonis_multi_rna_dna_path


### Now let's look at all of the PERMANOVA results together
adonis_multi_tax
adonis_multi_dna_path
adonis_multi_rna_path
adonis_multi_rna_dna_path


