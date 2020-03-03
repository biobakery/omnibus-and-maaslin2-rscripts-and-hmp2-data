# Maaslin2 per-feature testing using HMP2 week 0 subset data.

#Bring in the data needed. If you still have all of these loaded in your environment, then you can skip these steps.
#####
# MGX species
#####
# You can either download the file from bitbucket page or do it with the below code:
download.file("https://bitbucket.org/biobakery/omnibus-and-maaslin2-rscripts-and-hmp2-data/downloads/taxonomic_profiles_pcl_week0.csv", "./Data/taxonomic_profiles_pcl_week0.csv") # Download the mgx species data and put it into the data directory

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
# species[] = as.data.frame(sapply(species, function(x) as.numeric(as.character(x))))
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

#####
# MGX pathways
#####
# You can either download the file from bitbucket page or do it with the below code:
download.file("https://bitbucket.org/biobakery/omnibus-and-maaslin2-rscripts-and-hmp2-data/downloads/dna_pathabundance_relab_pcl_week0.csv", "./Data/dna_pathabundance_relab_pcl_week0.csv")

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

#####
# MTX pathways
#####
# You can either download the file from bitbucket page or do it with the below code:
download.file("https://bitbucket.org/biobakery/omnibus-and-maaslin2-rscripts-and-hmp2-data/downloads/rna_pathabundance_relab_pcl_week0.csv", "./Data/rna_pathabundance_relab_pcl_week0.csv")

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


#####
# RNA/DNA pathway ratios
#####
# You can either download the file from bitbucket page or do it with the below code:
download.file("https://bitbucket.org/biobakery/omnibus-and-maaslin2-rscripts-and-hmp2-data/downloads/rna_dna_path_relative_expression_week0.csv", "./Data/rna_dna_path_relative_expression_week0.csv")

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


###################################################################################################

### Feature-wise testing using MaAsLin2
# Load the package
library(Maaslin2)

#####
# MGX species

# Run Maaslin2
#defaults are tss normalization (turning that off since I already did that), log transform, lm test, max_significance=0.25,
#min_abundance=0, min_prevalence=0.1 (changing this to 0 because I already filtered). 
Maaslin2(species_filt, metadata, "maaslin2_tax", 
         normalization = "NONE", 
         min_abundance = 0, 
         min_prevalence = 0)

# Let's look over the results in the output folder
# Now that we've looked at it, we can see that diagnosis needs to be altered to assign nonIBD as the reference level.
# We can do this by adding an "a_" to the beginning in the metadata file
metadata$diagnosis = gsub("nonIBD", "a_nonIBD", metadata$diagnosis)
str(metadata) # looks good to go
# Now let's remove the old directory and rerun Maaslin2, we need to do this in terminal
# copy and paste this command in terminal within RStudio:

##### THIS IS TO BE RUN IN TERMINAL, NOT CONSOLE #####
rm -r maaslin2_tax

# Rerun Maaslin2
Maaslin2(species_filt, metadata, "maaslin2_tax", 
         normalization = "NONE", 
         min_abundance = 0, 
         min_prevalence = 0)

#Running all variables as univariate
dir.create(paste0(getwd(), '/maaslin2_tax_univ'))
for (i in 1:length(metadata))
{
  Maaslin2(species_filt, metadata[i], 
           output = paste0("maaslin2_tax_univ/", names(metadata)[i]), 
           normalization = "NONE", 
           min_abundance = 0, 
           min_prevalence = 0)
}


#####
# MGX pathways

# Run Maaslin2
#defaults are tss normalization (turning that off since I already did that), log transform, lm test, max_significance=0.25,
#min_abundance=0, min_prevalence=0.1 (changing this to 0 because I already filtered). 
Maaslin2(dna_path_unstratified_filt, metadata, "maaslin2_dna_path", 
         normalization = "NONE", 
         min_abundance = 0, 
         min_prevalence = 0)

#Running all variables as univariate
dir.create(paste0(getwd(), '/maaslin2_dna_path_univ'))
for (i in 1:length(metadata))
{
  Maaslin2(dna_path_unstratified_filt, metadata[i], 
           output = paste0("maaslin2_dna_path_univ/", names(metadata)[i]), 
           normalization = "NONE", 
           min_abundance = 0, 
           min_prevalence = 0)
}


#####
# MTX pathways

# Need to change diagnosis levels in metadata_rna_dna
metadata_rna$diagnosis = gsub("nonIBD", "a_nonIBD", metadata_rna$diagnosis)
str(metadata_rna) # looks good to go
# Run Maaslin2
#defaults are tss normalization (turning that off since I already did that), log transform, lm test, max_significance=0.25,
#min_abundance=0, min_prevalence=0.1 (changing this to 0 because I already filtered). 
Maaslin2(rna_path_unstratified_filt, metadata_rna, "maaslin2_rna_path", 
         normalization = "NONE", 
         min_abundance = 0, 
         min_prevalence = 0)

#Running all variables as univariate
dir.create(paste0(getwd(), '/maaslin2_rna_path_univ'))
for (i in 1:length(metadata_rna))
{
  Maaslin2(rna_path_unstratified_filt, metadata_rna[i], 
           output = paste0("maaslin2_rna_path_univ/", names(metadata_rna)[i]), 
           normalization = "NONE", 
           min_abundance = 0, 
           min_prevalence = 0)
}


#####
# MTX/MGX ratios pathways

# Need to change diagnosis levels in metadata_rna_dna
metadata_rna_dna$diagnosis = gsub("nonIBD", "a_nonIBD", metadata_rna_dna$diagnosis)
str(metadata_rna_dna) # looks good to go
# Run Maaslin2
#defaults are tss normalization (turning that off since I already did that), log transform, lm test, max_significance=0.25,
#min_abundance=0, min_prevalence=0.1 (changing this to 0 because I already filtered). 
# Changing transform to none, because we already transformed.
Maaslin2(rna_dna_path_unstratified_filt_log, metadata_rna_dna, "maaslin2_rna_dna_path", 
         transform = "NONE", 
         normalization = "NONE", 
         min_abundance = 0, 
         min_prevalence = 0)

#Running all variables as univariate
dir.create(paste0(getwd(), '/maaslin2_rna_dna_path_univ'))
for (i in 1:length(metadata_rna_dna))
{
  Maaslin2(rna_dna_path_unstratified_filt_log, metadata_rna_dna[i], 
           output = paste0("maaslin2_rna_dna_path_univ/", names(metadata_rna_dna)[i]), 
           transform = "NONE", 
           normalization = "NONE", 
           min_abundance = 0, 
           min_prevalence = 0)
}
