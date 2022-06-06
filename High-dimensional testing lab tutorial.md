# High-dimensional testing lab tutorial

- [Omnibus tests on MGX and MTX data from the HMP2 baseline only data](#omnibus-tests-on-mgx-and-mtx-data-from-the-hmp2-baseline-only-data)
  * [MGX taxonomy](#mgx-taxonomy)
    + [Feature table and metadata table creation and formatting](#feature-table-and-metadata-table-creation-and-formatting)
    + [Alpha diversity](#alpha-diversity)
      - [Univariable](#univariable)
      - [Multivariable](#multivariable)
    + [Beta Diversity](#beta-diversity)
      - [Univariable](#univariable)
      - [Multivariable](#multivariable)
      - [Pairwise Bray-Curtis comparisons](#pairwise-bray-curtis-comparisons)
  * [MGX pathways](#mgx-pathways)
    + [Feature table and metadata table creation and formatting](#feature-table-and-metadata-table-creation-and-formatting-1)
    + [Beta Diversity PERMANOVA tests using Bray-Curtis Dissimilarity](#beta-diversity-permanova-tests-using-bray-curtis-dissimilarity)
      - [Univariable](#univariable-2)
      - [Multivariable](#multivariable-2)
      - [Pairwise Bray-Curtis comparisons](#pairwise-bray-curtis-comparisons-1)
  * [MTX pathways](#mtx-pathways)
    + [Feature table and metadata table creation and formatting](#feature-table-and-metadata-table-creation-and-formatting-2)
    + [Beta Diversity PERMANOVA tests using Bray-Curtis Dissimilarity](#beta-diversity-permanova-tests-using-bray-curtis-dissimilarity-1)
      - [Univariable](#univariable-3)
      - [Multivariable](#multivariable-3)
      - [Pairwise Bray-Curtis comparisons](#pairwise-bray-curtis-comparisons-2)
  * [RNA DNA pathway ratios](#rna-dna-pathway-ratios)
    + [Feature table and metadata table creation and formatting](#feature-table-and-metadata-table-creation-and-formatting-3)
    + [Beta Diversity PERMANOVA tests using Euclidean distances](#beta-diversity-permanova-tests-using-euclidean-distances)
      - [Univariable](#univariable-4)
      - [Multivariable](#multivariable-4)
      - [Pairwise Euclidean comparisons](#pairwise-euclidean-comparisons)
- [Mantel tests](#mantel-tests)

# Omnibus tests on MGX and MTX data from the HMP2 baseline only data

Set the working directory:
```
setwd("Tutorials/highdimtesting")
```

Load the R packages needed:
```
library(vegan)
```

* Question: What are omnibus tests and how do they differ from featurewise tests? How might these complement each other?

## MGX taxonomy

### Feature table and metadata table creation and formatting

Download the MGX taxonomy relative abundance data and put it into the data directory:
```
download.file("https://raw.githubusercontent.com/biobakery/omnibus-and-maaslin2-rscripts-and-hmp2-data/master/taxonomic_profiles_pcl_week0.csv", "./Data/taxonomic_profiles_pcl_week0.csv")
```

Read the taxonomic relative abundance data into your R environment:
```
tax = read.csv(file = "./Data/taxonomic_profiles_pcl_week0.csv", header = T, row.names = 1, check.names = FALSE)
```

Take a look at `tax`:
```
tax[1:10,1:10]
```

* Question: Why not use head(tax)? Try it if you are unsure.

Check the dimensions:
```
dim(tax)
```
```
[1]   96 1484
```

Extract the metadata from `tax`, the first 5 columns in this file:
```
metadata = tax[1:5]
```

* Question: What are the advantages and disadvantages of storing and accessing data like this?  

Take a look at `metadata`:
```
head(metadata)
```
```
            site_name    sex  race consent_age diagnosis
CSM67UH7 Cedars-Sinai   Male White          69    nonIBD
CSM79HHW Cedars-Sinai Female White          36        CD
HSM67VDT   Cincinnati   Male White          11    nonIBD
HSM6XRQB   Cincinnati Female White          16        CD
HSM6XRR3   Cincinnati   Male White          13        CD
HSM7J4LP   Cincinnati Female White          12        CD
```

Check the structure:
```
str(metadata)
```
```
'data.frame':	96 obs. of  5 variables:
 $ site_name  : chr  "Cedars-Sinai" "Cedars-Sinai" "Cincinnati" "Cincinnati" ...
 $ sex        : chr  "Male" "Female" "Male" "Female" ...
 $ race       : chr  "White" "White" "White" "White" ...
 $ consent_age: int  69 36 11 16 13 12 16 14 30 57 ...
 $ diagnosis  : chr  "nonIBD" "CD" "nonIBD" "CD" ...
 ```

Check for NAs in `metadata`, since these will cause issues later in PERMANOVA tests:
```
sapply(metadata, function(x) sum(is.na(x)))
```
```
site_name         sex        race consent_age   diagnosis 
        0           0           0           6           0 
```

Age has 6 NA.

If this was a discrete variable we could classify the NAs as "unknown" and keep them in the model, but since age is a continuous variable typically we would either remove those from the data or impute. 

* Question: What is the main drawback of keeping NA values for discrete variables? Is there a case where this is totally justified?

In this case, let's impute with the median in order to not remove samples:
```
metadata$consent_age[is.na(metadata$consent_age)] = median(metadata$consent_age, na.rm = T)
```

Extract species data from `tax` and transpose the dataframe so that taxa are in rows: 
```
species = as.data.frame(t(tax[6:ncol(tax)]))
```
Note that transpose returns a matrix, so we need to coerce it back into a dataframe.

Check that `species` is all numeric values:
```
all(sapply(species, is.numeric))
```

Take a look at the output:
```
species[1:8,1:2]
```
```
                                                                                      CSM67UH7
k__Archaea                                                                                                                                                                                        0
k__Archaea|p__Euryarchaeota                                                                                                                                                                       0
k__Archaea|p__Euryarchaeota|c__Methanobacteria                                                                                                                                                    0
k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales                                                                                                                              0
k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae                                                                                                       0
k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae|g__Methanobrevibacter                                                                                 0
k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae|g__Methanobrevibacter|s__Methanobrevibacter_smithii                                                   0
k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae|g__Methanobrevibacter|s__Methanobrevibacter_smithii|t__Methanobrevibacter_smithii_unclassified        0
                                                                                                                                                                                           CSM79HHW
k__Archaea                                                                                                                                                                                        0
k__Archaea|p__Euryarchaeota                                                                                                                                                                       0
k__Archaea|p__Euryarchaeota|c__Methanobacteria                                                                                                                                                    0
k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales                                                                                                                              0
k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae                                                                                                       0
k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae|g__Methanobrevibacter                                                                                 0
k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae|g__Methanobrevibacter|s__Methanobrevibacter_smithii                                                   0
k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae|g__Methanobrevibacter|s__Methanobrevibacter_smithii|t__Methanobrevibacter_smithii_unclassified        0
```

As we can see here, `species` currently has all taxonomic levels for each microbe. For now, we want to only keep the species level information.

First, grep the rows that do not include strain stratifications and store them in a new temporary vector:
```
tmp_ind = grep("\\|t__", rownames(species), invert = T)
```

Create a new dataframe with only those row numbers in `tmp_ind`:
```
tmp = species[tmp_ind,]
```

Verify that strains were removed by peeking at the output:
```
row.names(tmp)[1:10]
```
```
  [1] "k__Archaea"                                                                                                                                                                           
  [2] "k__Archaea|p__Euryarchaeota"                                                                                                                                                          
  [3] "k__Archaea|p__Euryarchaeota|c__Methanobacteria"                                                                                                                                       
  [4] "k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales"                                                                                                                 
  [5] "k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae"                                                                                          
  [6] "k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae|g__Methanobrevibacter"                                                                    
  [7] "k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae|g__Methanobrevibacter|s__Methanobrevibacter_smithii"                                      
  [8] "k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae|g__Methanobrevibacter|s__Methanobrevibacter_unclassified"                                 
  [9] "k__Bacteria"                                                                                                                                                                          
 [10] "k__Bacteria|p__Acidobacteria"
```

* Question: What will happen if there is no strain level information in the starting file?

Now, let's select only the species level stratifications, removing all taxonomic levels before it:
```
tmp_ind = grep("\\|s__", rownames(tmp))
```

Create a new dataframe with only those row numbers in `tmp_ind`. Make sure to select from `tmp`, not `species`:
```
species = tmp[tmp_ind,]
```

Check the output to make sure that we only have species level stratifications:
```
row.names(species)[1:10]
```
```
  [1] "k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae|g__Methanobrevibacter|s__Methanobrevibacter_smithii"                                      
  [2] "k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae|g__Methanobrevibacter|s__Methanobrevibacter_unclassified"                                 
  [3] "k__Bacteria|p__Acidobacteria|c__Acidobacteriia|o__Acidobacteriales|f__Acidobacteriaceae|g__Granulicella|s__Granulicella_unclassified"                                                 
  [4] "k__Bacteria|p__Acidobacteria|c__Acidobacteriia|o__Acidobacteriales|f__Acidobacteriaceae|g__Terriglobus|s__Terriglobus_unclassified"                                                   
  [5] "k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinobaculum|s__Actinobaculum_schaalii"                                                    
  [6] "k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinobaculum|s__Actinobaculum_unclassified"                                                
  [7] "k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces|s__Actinomyces_europaeus"                                                       
  [8] "k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces|s__Actinomyces_graevenitzii"                                                    
  [9] "k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces|s__Actinomyces_odontolyticus"                                                   
 [10] "k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces|s__Actinomyces_turicensis"    
```

Remove temp files to clear up space:
```
rm(tmp, tmp_ind)
```
This isn't usually necessary, but R stores objects in RAM (virtual memory), so with very large datasets or small computing environments, it might matter.

Now we have only species level information in `species`, but let's make the names shorter for display by trimming off all other taxonomic information:
```
row.names(species) = gsub(".*\\|", "", row.names(species))
```

See if the row names look as expected:
```
row.names(species)[1:6]
```
```
  [1] "s__Methanobrevibacter_smithii"                 "s__Methanobrevibacter_unclassified"            "s__Granulicella_unclassified"                 
  [4] "s__Terriglobus_unclassified"                   "s__Actinobaculum_schaalii"                     "s__Actinobaculum_unclassified"                
```

Now that we have successfully truncated the names in `species`, let's check the sample sums (colSums) to make sure they are in proportion format (0-1) and are all ~1:
```
colSums(species)
```
```
   CSM67UH7    CSM79HHW    HSM67VDT    HSM6XRQB    HSM6XRR3    HSM7J4LP 
  0.9999752   1.0000004   1.0000005   0.9993355   0.9996256   1.0000003 
...
```

It is often useful to filter out low prevalence/abundance features in 'omics data, so we will make a filtered copy of the species table.

Check the dimensions of species pre-filtering:
```
dim(species)
```
```
[1] 572  96
```

Filter:
```
species_filt = species[apply(species, 1, function(x) sum(x > 0.0001) > 0.1 * ncol(species)), ]
```

* Question: What exactly is this command doing? What does 0.0001 represent? 0.1?
* Note that you can break the command apart to poke at the output, e.g. run `apply(species, 1, function(x) sum(x > 0.0001)` first.

Check the dimensions of species post-filtering:
```
dim(species_filt)
```
```
[1] 115  96
```

Let's transpose the dataframes for easier use downstream, making the rows be the samples, like in the metadata:
```
species_filt = as.data.frame(t(species_filt), check.names = F)
species = as.data.frame(t(species), check.names = F)
```
Note that `check.names = F` is supplied so that invalid or repeated column names are not flagged as errors. 
This is useful when you want to make sure your output always matches the original input file formating, but can be dangerous.  
The `make.names()` function can be used to sanitize character vectors.

Check to make sure the sample names are in the same order in both the metadata and the species dataframes:
```
all(row.names(metadata) == row.names(species_filt))
all(row.names(metadata) == row.names(species))
```

* Question: Why perform this type of filtering? Are there situations were it might not be appropriate?

>NOTE: These formatted files are also located in the Tutorials/highdimtesting directory of the bioBakery image. To work with them from there assign them in R with the following code:
>```
>metadata = read.csv(file = "metadata.csv", header = T, row.names = 1, check.names = F)
>species = read.csv(file = "species.csv", header = T, row.names = 1, check.names = F)
>species_filt = read.csv(file = "species_filt.csv", header = T, row.names = 1, check.names = F)
>```

### Alpha diversity

* Question: What is alpha diversity, and why is it important? 

Using the unfiltered species data, calculate a common alpha diversity index:
```
alpha_div = as.data.frame(diversity(species, index = "shannon"))
```
You can check which indices the vegan package supports with `?diversity`. We often also use Inverse Simpson, which gives similar results. 
However, there are several more you might see in the literature, with varying pros and cons, e.g. Chao1.

Check the output:
```
head(alpha_div)
```
```
BLOOP
```

Rename the columns to shorter names:
```
names(alpha_div) = c("shannon")
```

Collate with the metadata, so that a single dataframe can be provided to upcoming statistical tests:
```
alpha_meta <- merge(alpha_div, metadata, by = "row.names")
```
Note that the `merge` function makes the old row names into new column. This might not be what you expect, but you can confirm with `?merge`.

Restore the row names (sample IDs) and get rid of the column generated by the merge:
```
row.names(alpha_meta) = alpha_meta$Row.names
alpha_meta$Row.names = NULL
```

Confirm this looks as expected:
```
head(alpha_meta)
```
```
BLOOP
```

#### Univariable

Alpha diversity is relatively well-behaved and generates a single value per sample, allowing for a number of statistical tests.
Linear models are straightforward and can accommodate a large variety of experimental designs, so we will go with that.

Run a linear model on diagnosis and return as an ANOVA table:
```
alpha_shannon_univ = anova(lm(Shannon ~ diagnosis, data = alpha_meta)) 
```

Can do this with a for loop to quickly get results for each metadata variable:
```
for (col in names(metadata)){
  alpha_shannon_univ = anova(lm(as.formula(paste("Shannon ~ ", col)), data = alpha_meta)) 
  print(col)
  print(alpha_shannon_univ)
}
```
```
BLOOP
```

* Question: Everyone likely knows what the Pr(>F) value means (the p-value), but what does the Df value tell us?

#### Multivariable

* Question: What are univariable and multivariable tests and what benefit do both provide?

Try a model with the age and diagnosis:
```
shannon_multi = anova(lm(Shannon ~ consent_age + diagnosis, data = alpha_met_df))
shannon_multi
```
```
BLOOP
```

If all variables in the file are used, there's a nice shorthand in R:
```
shannon_multi = anova(lm(Shannon ~ ., data = alpha_meta))
shannon_multi
```
```
BLOOP
```

* Question: This can be expanded to include random effects with the `lme4` package. Is there anything in this demo data that might indicate a mixed model would be useful?

### Beta Diversity

* Question: What is beta diversity, and why is it important? 

Using the filtered species table, calculate Bray-Curtis dissimilarity:
```
bray_species = vegdist(species_filt, method = "bray")
```

* Question: What makes Bray-Curtis dissimilarity a good index to use for this type of data? Do you know of other options?

#### Univariable

Unlike alpha diversity, beta diversity is defined between samples, and therefore, many common statistical tests are not usable. 

* Question: Why bother with beta diversity then?

PERMANOVA is the most common omnibus beta diversity test and is implemented as `adonis2` by the vegan package.
Note that `adonis` is now officially deprecated and older code that uses this function may act strangely. 

Run permanova and print the results:
```
set.seed(123)
adonis_diagnosis_tax = adonis2(bray_species ~ diagnosis, data = metadata)
adonis_diagnosis_tax
```
```
BLOOP
```
As you can see, the results look quite similar to the more familiar ANOVA test.

* Question: Why is a seed set before calling `adonis`? Try running it again without resetting the seed.

Can do this with a for loop:
```
for (col in names(metadata)){
  set.seed(123)
  adonis_univ <- adonis2(as.formula(paste("bray_species ~ ", col)), data = metadata)
  print(adonis_univ)
}
```
```
BLOOP
```

* Question: What does the R2 value tell us?

#### Multivariable

BLOOP:
```
adonis_multi_tax = adonis(bray_species ~ site_name + sex + race + consent_age + diagnosis, data = metadata)
adonis_multi_tax
```
```
BLOOP
```

Can use the same shorthand as with linear models:
```
adonis_multi_tax = adonis(bray_species ~ ., data = metadata)
adonis_multi_tax
```
```
BLOOP
```

## MGX pathways

### Feature table and metadata table creation and formatting

Download the mgx pathway relative abundance data and put it into the data directory:
```
download.file("https://raw.githubusercontent.com/biobakery/omnibus-and-maaslin2-rscripts-and-hmp2-data/master/dna_pathabundance_relab_pcl_week0.csv", "./Data/dna_pathabundance_relab_pcl_week0.csv")
```

Read the dna pathway relative abundance data into your R environment:
```
dna_path = read.csv(file = "./Data/dna_pathabundance_relab_pcl_week0.csv", header = T, row.names = 1, check.names = FALSE)
```

Check out `dna_path`:

Check the dimensions:
```
dim(dna_path)
```
```
[1]    96 10889
```

Remove metadata and keep only pathways and transpose the data:
```
dna_path = data.frame(t(dna_path[6:ncol(dna_path)]))
```

Check out `dna_path`:

Check out the structure:
```
str(dna_path)
```
```
'data.frame':	10884 obs. of  96 variables:
 $ CSM67UH7   : num  1.37e-02 7.18e-03 0.00 1.53e-04 1.62e-05 ...
 $ CSM79HHW   : num  0.017841 0 0 0.000633 0 ...
 $ HSM67VDT   : num  0.021241 0.000319 0 0.000311 0 ...
 $ HSM6XRQB   : num  0.01168 0 0 0.000518 0 ...
 $ HSM6XRR3   : num  0.014899 0.000186 0 0 0 ...
 $ HSM7J4LP   : num  0.00248 0 0 0 0 ...
 ...
```
Everything is numeric; good to go.

Check the output:
```
dna_path[1:8,1:2]
```
```
                                                            CSM67UH7    CSM79HHW
1CMET2-PWY: N10-formyl-tetrahydrofolate biosynthesis                                                0.013653600 0.017841100
1CMET2-PWY: N10-formyl-tetrahydrofolate biosynthesis|g__Akkermansia.s__Akkermansia_muciniphila      0.007178950 0.000000000
1CMET2-PWY: N10-formyl-tetrahydrofolate biosynthesis|g__Bacteroides.s__Bacteroides_barnesiae        0.000000000 0.000000000
1CMET2-PWY: N10-formyl-tetrahydrofolate biosynthesis|g__Bacteroides.s__Bacteroides_caccae           0.000152978 0.000632591
1CMET2-PWY: N10-formyl-tetrahydrofolate biosynthesis|g__Bacteroides.s__Bacteroides_cellulosilyticus 0.000016200 0.000000000
1CMET2-PWY: N10-formyl-tetrahydrofolate biosynthesis|g__Bacteroides.s__Bacteroides_clarus           0.000000000 0.000000000
1CMET2-PWY: N10-formyl-tetrahydrofolate biosynthesis|g__Bacteroides.s__Bacteroides_coprocola        0.000000000 0.000000000
1CMET2-PWY: N10-formyl-tetrahydrofolate biosynthesis|g__Bacteroides.s__Bacteroides_faecis           0.000000000 0.000000000
```

As we can see here, `dna_path` currently has the microbial species that contribute to each dna pathway, so now we need to only keep the dna pathway information.

grep the rows that do not include species stratifications and store them in a new temporary vector:
```
tmp.ind = grep("\\|.*", rownames(dna_path), invert = T)
```

Create a new dataframe with only those row numbers in `tmp.ind`:
```
dna_path_unstratified = dna_path[tmp.ind,]
```
Check the output:
```
row.names(dna_path_unstratified)[1:10]
```
```
  [1] "1CMET2-PWY: N10-formyl-tetrahydrofolate biosynthesis"                                                          
  [2] "3-HYDROXYPHENYLACETATE-DEGRADATION-PWY: 4-hydroxyphenylacetate degradation"                                    
  [3] "7ALPHADEHYDROX-PWY: cholate degradation (bacteria, anaerobic)"                                                 
  [4] "AEROBACTINSYN-PWY: aerobactin biosynthesis"                                                                    
  [5] "ALL-CHORISMATE-PWY: superpathway of chorismate metabolism"                                                     
  [6] "ALLANTOINDEG-PWY: superpathway of allantoin degradation in yeast"                                              
  [7] "ANAEROFRUCAT-PWY: homolactic fermentation"                                                                     
  [8] "ANAGLYCOLYSIS-PWY: glycolysis III (from glucose)"                                                              
  [9] "ARG+POLYAMINE-SYN: superpathway of arginine and polyamine biosynthesis"                                        
 [10] "ARGDEG-PWY: superpathway of L-arginine, putrescine, and 4-aminobutanoate degradation"
```

Remove temp files to clear up space:
```
rm(tmp.ind)
```

Now that we have sucessfully unstratified `dna_path`, let's check the sample sums (colSums) to make sure they are in proportion format (0-1) and are all ~1:
```
colSums(dna_path_unstratified)
```
```
   CSM67UH7    CSM79HHW    HSM67VDT    HSM6XRQB    HSM6XRR3    HSM7J4LP 
  1.0000000   1.0000005   1.0000002   0.9999997   1.0000002   0.9999997 
...
```

Filter for beta diversity:

Check the dimensions of dna pathways pre-filtering:
```
dim(dna_path_unstratified)
```
```
[1] 466  96
```
Filter:
```
dna_path_unstratified_filt = dna_path_unstratified[apply(dna_path_unstratified, 1, function(x) sum(x > 0.0001) > 0.1 * ncol(dna_path_unstratified)), ]
```
Check the dimensions of dna pathways post-filtering:
```
dim(dna_path_unstratified_filt)
```
```
[1] 309  96
```

Let's transpose the dataframes for easier use downstream, making the rows be the samples just like in the metadata:
```
dna_path_unstratified_filt = data.frame(t(dna_path_unstratified_filt), check.names = F)
dna_path_unstratified = data.frame(t(dna_path_unstratified), check.names = F)
```

Check to make sure the sample names are in the same order in both the metadata and the DNA pathways dataframe:
```
row.names(metadata) == row.names(dna_path_unstratified_filt)
```
```
 [1]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
[13]  TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE
[25]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
[37]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
[49]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
[61]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
[73]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
[85]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
```
There are two samples that do not match. Let's match them.

```
dna_path_unstratified_filt = dna_path_unstratified_filt[match(row.names(metadata), row.names(dna_path_unstratified_filt)),]
```
```
row.names(metadata) == row.names(dna_path_unstratified_filt)
```
```
 [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
[16] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
[31] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
[46] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
[61] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
[76] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
[91] TRUE TRUE TRUE TRUE TRUE TRUE
```
Good to go now.

>NOTE: These formatted files are also located in the Tutorials/highdimtesting directory of the bioBakery image. To work with them from there just assign them in R with the following code:
>```
>dna_path_unstratified_filt = read.csv(file = "dna_path_unstratified_filt.csv", header = T, row.names = 1, check.names = FALSE)
>```


### Beta Diversity PERMANOVA tests using Bray-Curtis Dissimilarity

Calculate Bray-Curtis dissimilarity:
```
bray_dna_path_unstratified = vegdist(dna_path_unstratified_filt, method = "bray")
```

#### Univariable

Can do this with a for loop:
```
for (col in names(metadata)){
  adonis_univ <- adonis(as.formula(paste("bray_dna_path_unstratified ~ ", col)), data = metadata)
  print(adonis_univ)
}
```
```
Call:
adonis(formula = as.formula(paste("bray_dna_path_unstratified ~ ",      col)), data = metadata) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
site_name  4    0.1763 0.044076  1.1084 0.04646  0.294
Residuals 91    3.6187 0.039765         0.95354       
Total     95    3.7950                  1.00000       

Call:
adonis(formula = as.formula(paste("bray_dna_path_unstratified ~ ",      col)), data = metadata) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs  MeanSqs F.Model     R2 Pr(>F)
sex        1    0.0383 0.038313 0.95867 0.0101  0.435
Residuals 94    3.7566 0.039964         0.9899       
Total     95    3.7950                  1.0000       

Call:
adonis(formula = as.formula(paste("bray_dna_path_unstratified ~ ",      col)), data = metadata) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
race       3    0.1211 0.040381  1.0112 0.03192  0.414
Residuals 92    3.6738 0.039933         0.96808       
Total     95    3.7950                  1.00000       

Call:
adonis(formula = as.formula(paste("bray_dna_path_unstratified ~ ",      col)), data = metadata) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
consent_age  1    0.0367 0.036663   0.917 0.00966  0.457
Residuals   94    3.7583 0.039982         0.99034       
Total       95    3.7950                  1.00000       

Call:
adonis(formula = as.formula(paste("bray_dna_path_unstratified ~ ",      col)), data = metadata) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
diagnosis  2    0.1327 0.066335  1.6845 0.03496  0.061 .
Residuals 93    3.6623 0.039379         0.96504         
Total     95    3.7950                  1.00000         
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

Alternatively, can do it one-by-one like above in the taxonomy section.


#### Multivariable

Can do it without being verbose:
```
adonis_multi_dna_path = adonis(bray_dna_path_unstratified ~ ., data = metadata)
adonis_multi_dna_path
```
```
Call:
adonis(formula = bray_dna_path_unstratified ~ ., data = metadata) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
site_name    4    0.1763 0.044076  1.1105 0.04646  0.273
sex          1    0.0311 0.031101  0.7836 0.00820  0.616
race         3    0.1243 0.041425  1.0437 0.03275  0.394
consent_age  1    0.0293 0.029283  0.7378 0.00772  0.620
diagnosis    2    0.1001 0.050034  1.2606 0.02637  0.222
Residuals   84    3.3339 0.039690         0.87852       
Total       95    3.7950                  1.00000       
```

* Question: Wait a minute! Diagnosis was near significant in the univariable model, but is now quite non-significant in the multivariable model. What happened?

Alternatively, can write out each variable in the model like above in the taxonomy section.

* Exercise: Run the same model a few times and compare the results.
    * Question: Are you seeing the same results across runs? If not, why is this?

#### Pairwise Bray-Curtis comparisons
Pairwise comparisons of diagnosis in a univariable model

```
bray_dna_path_unstratified_mat = as.matrix(bray_dna_path_unstratified)
pairwise.adonis(bray_dna_path_unstratified_mat, factors = metadata$diagnosis, p.adjust.m = "BH")
```
```
         pairs Df  SumsOfSqs  F.Model         R2 p.value p.adjusted sig
1 nonIBD vs CD  1 0.01580812 1.021681 0.01459093   0.343     0.3430    
2 nonIBD vs UC  1 0.05216168 3.608071 0.07129439   0.032     0.0960    
3     CD vs UC  1 0.04085920 2.340545 0.03235454   0.085     0.1275    
```


## MTX pathways

### Feature table and metadata table creation and formatting

Download the mtx pathway relative abundance data and put it into the data directory:
```
download.file("https://raw.githubusercontent.com/biobakery/omnibus-and-maaslin2-rscripts-and-hmp2-data/master/rna_pathabundance_relab_pcl_week0.csv", "./Data/rna_pathabundance_relab_pcl_week0.csv")
```

Read the rna pathway data into R environment:
```
rna_path = read.csv(file = "./Data/rna_pathabundance_relab_pcl_week0.csv", header = T, row.names = 1, check.names = FALSE)
```

Check out `rna_path`:

Check the dimensions:
```
dim(rna_path)
```
```
[1]   28 6066
```

* Question: What is a possible reason for this rna data having less samples than the dna data?

Remove metadata and keep only pathways and transpose the data:
```
rna_path = data.frame(t(rna_path[6:ncol(rna_path)]))
```

Check out `rna_path`:

Check out the structure:
```
str(rna_path)
```
```
'data.frame':	6061 obs. of  28 variables:
 $ HSM6XRQB   : num  0.0167 0 0 0 0 ...
 $ HSM7J4LP   : num  0 0 0 0 0 0 0 0 0 0 ...
 $ MSM5LLDI   : num  0.0198 0 0 0 0 ...
 $ MSM6J2Q1   : num  0.0133 0 0 0 0 ...
 $ MSM79H8D   : num  0.0186 0 0 0 0 ...
 $ MSM9VZL5   : num  0.015814 0 0 0.000274 0 ...
 $ MSM9VZMM   : num  0.0176 0 0 0 0 ...
 $ PSM6XBRK   : num  0.0149 0 0 0 0 ...
 $ PSM6XBSE   : num  0.00742 0 0 0 0 ...
 $ PSM7J4EF   : num  0.00555 0 0 0 0 ...
 $ PSMA265X   : num  0.0113 0 0 0 0 ...
 $ PSM6XBRK_TR: num  0.014 0 0 0 0 ...
 $ MSM5LLHR_P : num  0.00566 0 0 0 0 ...
 $ ESM5MEDZ_P : num  0.0109 0 0 0 0 ...
 $ MSM5LLIQ_P : num  0.0105 0 0 0 0 ...
 $ MSM5LLIC_P : num  0.01011 0 0 0.000189 0 ...
 $ MSM5LLHV_P : num  0.0178 0 0 0 0 ...
 $ CSM5FZ3T_P : num  0.0316 0 0 0 0 ...
 $ MSM5LLFG_P : num  0.0173 0 0 0 0 ...
 $ MSM5LLIS_P : num  0.00486 0 0 0 0 ...
 $ ESM5GEYY_P : num  0.0121 0 0 0 0 ...
 $ CSM5MCTZ_P : num  0.012069 0 0 0.000355 0 ...
 $ CSM5MCV1_P : num  0.00664 0 0 0.000394 0 ...
 $ HSM5MD87_P : num  0.009551 0 0 0.000383 0 ...
 $ CSM5MCU4_P : num  0.00802 0 0 0 0 ...
 $ HSM5MD82_P : num  0 0 0 0 0 0 0 0 0 0 ...
 $ HSM5MD8A_P : num  0.00209 0 0 0 0 ...
 $ CSM6J2H9_P : num  0.0271 0 0 0 0 ...
```
Everything is numeric; good to go.

Check the output:
```
rna_path[1:8,1:2]
```
```
                                                            HSM6XRQB HSM7J4LP
1CMET2-PWY: N10-formyl-tetrahydrofolate biosynthesis                                                0.0166822        0
1CMET2-PWY: N10-formyl-tetrahydrofolate biosynthesis|g__Akkermansia.s__Akkermansia_muciniphila      0.0000000        0
1CMET2-PWY: N10-formyl-tetrahydrofolate biosynthesis|g__Bacteroides.s__Bacteroides_barnesiae        0.0000000        0
1CMET2-PWY: N10-formyl-tetrahydrofolate biosynthesis|g__Bacteroides.s__Bacteroides_caccae           0.0000000        0
1CMET2-PWY: N10-formyl-tetrahydrofolate biosynthesis|g__Bacteroides.s__Bacteroides_cellulosilyticus 0.0000000        0
1CMET2-PWY: N10-formyl-tetrahydrofolate biosynthesis|g__Bacteroides.s__Bacteroides_clarus           0.0000000        0
1CMET2-PWY: N10-formyl-tetrahydrofolate biosynthesis|g__Bacteroides.s__Bacteroides_coprocola        0.0000000        0
1CMET2-PWY: N10-formyl-tetrahydrofolate biosynthesis|g__Bacteroides.s__Bacteroides_faecis           0.0000000        0
```

Minimize the metadata to just the samples available in these data:

Check the dimensions of `metadata` pre-subsetting:
```
dim(metadata)
```
```
[1] 96  5
```

Subset metadata and assign it to a new dataframe:
```
metadata_rna = subset(metadata, row.names(metadata) %in% names(rna_path))
```

Check the dimensions of the subsetted metadata:
```
dim(metadata_rna)
```
```
[1] 28  5
```

Check the output:
```
metadata_rna
```
```
                 site_name    sex                      race consent_age diagnosis
HSM6XRQB        Cincinnati Female                     White          16        CD
HSM7J4LP        Cincinnati Female                     White          12        CD
MSM5LLDI               MGH Female                     White          30        CD
MSM6J2Q1               MGH   Male                     White          28    nonIBD
MSM79H8D               MGH   Male                     White          61    nonIBD
MSM9VZL5               MGH   Male                     White          21    nonIBD
MSM9VZMM               MGH Female                     White          23    nonIBD
PSM6XBRK    MGH Pediatrics   Male                     White          16        CD
PSM6XBRK_TR MGH Pediatrics   Male                     White          16        CD
PSM6XBSE    MGH Pediatrics Female                     White          16        UC
PSM7J4EF    MGH Pediatrics   Male                     White          15        CD
PSMA265X    MGH Pediatrics Female                     White          16        UC
MSM5LLHR_P             MGH   Male                     White          26        CD
ESM5MEDZ_P           Emory Female        More than one race           8        CD
MSM5LLIQ_P             MGH Female                     White          21        UC
MSM5LLIC_P             MGH Female                     White          43        CD
MSM5LLHV_P             MGH Female                     White          38        UC
CSM5FZ3T_P    Cedars-Sinai Female                     White          76        CD
MSM5LLFG_P             MGH Female                     White          24        CD
MSM5LLIS_P             MGH   Male                     Other          41        CD
ESM5GEYY_P           Emory   Male Black or African American          20        CD
CSM5MCTZ_P    Cedars-Sinai   Male                     White          32        UC
CSM5MCV1_P    Cedars-Sinai Female                     White          20        CD
HSM5MD87_P      Cincinnati   Male Black or African American          13        UC
CSM5MCU4_P    Cedars-Sinai Female                     White          53        CD
HSM5MD82_P      Cincinnati Female                     White           6    nonIBD
HSM5MD8A_P      Cincinnati Female                     White          13    nonIBD
CSM6J2H9_P             MGH Female                     White          40    nonIBD
```

As we can see above, `rna_path` currently has the microbial species that contribute to each rna pathway, so now we need to only keep the rna pathway information.

grep the rows that do not include species stratifications and store them in a new temporary vector:
```
tmp.ind = grep("\\|.*", rownames(rna_path), invert = T)
```

Create a new dataframe with only those row numbers in `tmp.ind`:
```
rna_path_unstratified = rna_path[tmp.ind,]
```
Check the output:
```
row.names(rna_path_unstratified)[1:10]
```
```
  [1] "1CMET2-PWY: N10-formyl-tetrahydrofolate biosynthesis"                                                          
  [2] "3-HYDROXYPHENYLACETATE-DEGRADATION-PWY: 4-hydroxyphenylacetate degradation"                                    
  [3] "7ALPHADEHYDROX-PWY: cholate degradation (bacteria, anaerobic)"                                                 
  [4] "AEROBACTINSYN-PWY: aerobactin biosynthesis"                                                                    
  [5] "ALL-CHORISMATE-PWY: superpathway of chorismate metabolism"                                                     
  [6] "ALLANTOINDEG-PWY: superpathway of allantoin degradation in yeast"                                              
  [7] "ANAEROFRUCAT-PWY: homolactic fermentation"                                                                     
  [8] "ANAGLYCOLYSIS-PWY: glycolysis III (from glucose)"                                                              
  [9] "ARG+POLYAMINE-SYN: superpathway of arginine and polyamine biosynthesis"                                        
 [10] "ARGDEG-PWY: superpathway of L-arginine, putrescine, and 4-aminobutanoate degradation" 
```

Remove temp files to clear up space:
```
rm(tmp.ind)
```

Now that we have sucessfully unstratified `rna_path`, let's check the sample sums (colSums) to make sure they are in proportion format (0-1) and are all ~1:
```
colSums(rna_path_unstratified)
```
```
   HSM6XRQB    HSM7J4LP    MSM5LLDI    MSM6J2Q1    MSM79H8D    MSM9VZL5 
  1.0000001   0.9999992   0.9999997   1.0000000   1.0000001   0.9999998 
...
```

Filter for beta diversity:

Check the dimensions of rna pathways pre-filtering:
```
dim(rna_path_unstratified)
```
```
[1] 430  28
```
Filter:
```
rna_path_unstratified_filt = rna_path_unstratified[apply(rna_path_unstratified, 1, function(x) sum(x > 0.0001) > 0.1 * ncol(rna_path_unstratified)), ]
```
Check the dimensions of rna pathways post-filtering:
```
dim(rna_path_unstratified_filt)
```
```
[1] 234  28
```

Let's transpose the dataframes for easier use downstream, making the rows be the samples just like in the metadata:
```
rna_path_unstratified_filt = data.frame(t(rna_path_unstratified_filt), check.names = F)
rna_path_unstratified = data.frame(t(rna_path_unstratified), check.names = F)
```

Check to make sure the sample names are in the same order in both the metadata and the RNA pathways dataframe:
```
row.names(metadata_rna) == row.names(rna_path_unstratified_filt)
```
```
 [1]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE
[13]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
[25]  TRUE  TRUE  TRUE  TRUE
```
They do not match. Let's match them.

```
rna_path_unstratified_filt = rna_path_unstratified_filt[match(row.names(metadata_rna), row.names(rna_path_unstratified_filt)),]
```
```
row.names(metadata_rna) == row.names(rna_path_unstratified_filt)
```
```
 [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
[16] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
```
Good to go now.

>NOTE: These formatted files are also located in the Tutorials/highdimtesting directory of the bioBakery image. To work with them from there just assign them in R with the following code:
>```
>metadata_rna = read.csv(file = "metadata_rna.csv", header = T, row.names = 1, check.names = FALSE)
>rna_path_unstratified_filt = read.csv(file = "rna_path_unstratified_filt.csv", header = T, row.names = 1, check.names = FALSE)
>```

### Beta Diversity PERMANOVA tests using Bray-Curtis Dissimilarity

Calculate Bray-Curtis dissimilarity:
```
bray_rna_path_unstratified = vegdist(rna_path_unstratified_filt, method = "bray")
```

#### Univariable

Can do this with a for loop:
```
for (col in names(metadata)){
  adonis_univ <- adonis(as.formula(paste("bray_rna_path_unstratified ~ ", col)), data = metadata)
  print(adonis_univ)
}
```
```
Error in G * t(hat) : non-conformable arrays
```


* Question: Why does this result in an error? 
    * Hint: Check the formula call and remember what we subsetted for this rna data.


Alternatively, can do it one-by-one like in the taxonomy section.


Solution to the for loop question:
```
for (col in names(metadata_rna)){
  adonis_univ <- adonis(as.formula(paste("bray_rna_path_unstratified ~ ", col)), data = metadata_rna)
  print(adonis_univ)
}
```
```
Call:
adonis(formula = as.formula(paste("bray_rna_path_unstratified ~ ",      col)), data = metadata_rna) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
site_name  4   0.33933 0.084834  1.0482 0.15419  0.369
Residuals 23   1.86136 0.080929         0.84581       
Total     27   2.20069                  1.00000       

Call:
adonis(formula = as.formula(paste("bray_rna_path_unstratified ~ ",      col)), data = metadata_rna) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs  MeanSqs F.Model     R2 Pr(>F)
sex        1    0.0823 0.082297  1.0101 0.0374  0.435
Residuals 26    2.1184 0.081477         0.9626       
Total     27    2.2007                  1.0000       

Call:
adonis(formula = as.formula(paste("bray_rna_path_unstratified ~ ",      col)), data = metadata_rna) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
race       3   0.13184 0.043945 0.50979 0.05991  0.909
Residuals 24   2.06886 0.086202         0.94009       
Total     27   2.20069                  1.00000       

Call:
adonis(formula = as.formula(paste("bray_rna_path_unstratified ~ ",      col)), data = metadata_rna) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
consent_age  1   0.08368 0.083678  1.0277 0.03802  0.318
Residuals   26   2.11701 0.081424         0.96198       
Total       27   2.20069                  1.00000       

Call:
adonis(formula = as.formula(paste("bray_rna_path_unstratified ~ ",      col)), data = metadata_rna) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
diagnosis  2   0.11062 0.055308 0.66155 0.05026  0.921
Residuals 25   2.09008 0.083603         0.94974       
Total     27   2.20069                  1.00000        
```

#### Multivariable

Can do it without being verbose:
```
adonis_multi_rna_path = adonis(bray_rna_path_unstratified ~ ., data = metadata_rna)
adonis_multi_rna_path
```
```
Call:
adonis(formula = bray_rna_path_unstratified ~ ., data = metadata_rna) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
site_name    4   0.33933 0.084834 0.94530 0.15419  0.489
sex          1   0.05461 0.054609 0.60851 0.02481  0.888
race         3   0.16446 0.054820 0.61086 0.07473  0.819
consent_age  1   0.06739 0.067385 0.75088 0.03062  0.558
diagnosis    2   0.13903 0.069515 0.77461 0.06318  0.811
Residuals   16   1.43587 0.089742         0.65246       
Total       27   2.20069                  1.00000       
```

Alternatively, can write out each variable in the model like in the taxonomy section.

#### Pairwise Bray-Curtis comparisons
Pairwise comparisons of diagnosis in a univariable model

```
bray_rna_path_unstratified_mat = as.matrix(bray_rna_path_unstratified)
pairwise.adonis(bray_rna_path_unstratified_mat, factors = metadata_rna$diagnosis, p.adjust.m = "BH")
```
```
         pairs Df   SumsOfSqs   F.Model         R2 p.value p.adjusted sig
1 CD vs nonIBD  1 0.018338280 1.0814873 0.05130033   0.311      0.804    
2     CD vs UC  1 0.011695096 0.6518523 0.03317002   0.536      0.804    
3 nonIBD vs UC  1 0.003735382 0.3904319 0.03427718   0.822      0.822    
```


## RNA DNA pathway ratios

### Feature table and metadata table creation and formatting

Download the RNA/DNA pathway ratio data and put it into the data directory:
```
download.file("https://raw.githubusercontent.com/biobakery/omnibus-and-maaslin2-rscripts-and-hmp2-data/master/rna_dna_path_relative_expression_week0.csv", "./Data/rna_dna_path_relative_expression_week0.csv")
```

Read the RNA/DNA pathway data into R environment:
```
rna_dna_path = read.csv(file = "./Data/rna_dna_path_relative_expression_week0.csv", header = T, row.names = 1, check.names = FALSE)
```
Check out `rna_dna_path`:

Check the dimensions:
```
dim(rna_dna_path)
```
```
[1]    12 10598
```

* Question: Why might this RNA/DNA pathway ratios dataset have less samples than the previous RNA pathway relative abundance dataset, 12 instead of 28?

This data file does not contain metadata like the others, so let's transpose the data:
```
rna_dna_path = data.frame(t(rna_dna_path))
```

Check out `rna_dna_path`:

Check out the structure:
```
str(rna_dna_path)
```
```
'data.frame':	10598 obs. of  12 variables:
 $ MSM5LLDI   : num  0.973 0 0 0 0 ...
 $ MSM9VZL5   : num  1.23 0 0 4.08 0 ...
 $ PSM6XBRK   : num  0.841 0 0 0 0 ...
 $ PSMA265X   : num  0.773 0 0 0 0 ...
 $ MSM79H8D   : num  1.08 0 0 0 0 ...
 $ PSM7J4EF   : num  0.697 0 0 0 0 ...
 $ MSM9VZMM   : num  1.06 0 0 0 0 ...
 $ HSM7J4LP   : num  0 0 0 0 0 0 0 0 0 0 ...
 $ PSM6XBSE   : num  0.822 0 0 0 0 ...
 $ MSM6J2Q1   : num  1.13 0 0 0 0 ...
 $ PSM6XBRK_TR: num  0.839 0 0 0 0 ...
 $ HSM6XRQB   : num  1.43 0 0 0 0 ...
```
Everything is numeric; good to go.

Check the output:
```
rna_dna_path[1:8,1:2]
```
```
                                                                                                     MSM5LLDI MSM9VZL5
1CMET2-PWY: N10-formyl-tetrahydrofolate biosynthesis                                                0.9733041 1.231584
1CMET2-PWY: N10-formyl-tetrahydrofolate biosynthesis|g__Akkermansia.s__Akkermansia_muciniphila      0.0000000 0.000000
1CMET2-PWY: N10-formyl-tetrahydrofolate biosynthesis|g__Bacteroides.s__Bacteroides_barnesiae        0.0000000 0.000000
1CMET2-PWY: N10-formyl-tetrahydrofolate biosynthesis|g__Bacteroides.s__Bacteroides_caccae           0.0000000 4.075788
1CMET2-PWY: N10-formyl-tetrahydrofolate biosynthesis|g__Bacteroides.s__Bacteroides_cellulosilyticus 0.0000000 0.000000
1CMET2-PWY: N10-formyl-tetrahydrofolate biosynthesis|g__Bacteroides.s__Bacteroides_clarus           0.0000000 0.000000
1CMET2-PWY: N10-formyl-tetrahydrofolate biosynthesis|g__Bacteroides.s__Bacteroides_coprocola        0.0000000 0.000000
1CMET2-PWY: N10-formyl-tetrahydrofolate biosynthesis|g__Bacteroides.s__Bacteroides_faecis           0.0000000 0.000000
```

Minimize the metadata to just the samples available in these data:

Check the dimensions of `metadata` pre-subsetting:
```
dim(metadata)
```
```
[1] 96  5
```

Subset metadata and assign it to a new dataframe:
```
metadata_rna_dna = subset(metadata, row.names(metadata) %in% names(rna_dna_path))
```

Check the dimensions of the subsetted metadata:
```
dim(metadata_rna_dna)
```
```
[1] 12  5
```

Check the output:
```
metadata_rna_dna
```
```
                 site_name    sex  race consent_age diagnosis
HSM6XRQB        Cincinnati Female White          16        CD
HSM7J4LP        Cincinnati Female White          12        CD
MSM5LLDI               MGH Female White          30        CD
MSM6J2Q1               MGH   Male White          28    nonIBD
MSM79H8D               MGH   Male White          61    nonIBD
MSM9VZL5               MGH   Male White          21    nonIBD
MSM9VZMM               MGH Female White          23    nonIBD
PSM6XBRK    MGH Pediatrics   Male White          16        CD
PSM6XBRK_TR MGH Pediatrics   Male White          16        CD
PSM6XBSE    MGH Pediatrics Female White          16        UC
PSM7J4EF    MGH Pediatrics   Male White          15        CD
PSMA265X    MGH Pediatrics Female White          16        UC
```

As we can see above, `rna_dna_path` currently has the microbial species that contribute to each rna pathway, so now we need to only keep the RNA/DNA pathway information.

grep the rows that do not include species stratifications and store them in a new temporary vector:
```
tmp.ind = grep("\\|.*", rownames(rna_dna_path), invert = T)
```

Create a new dataframe with only those row numbers in `tmp.ind`:
```
rna_dna_path_unstratified = rna_dna_path[tmp.ind,]
```
Check the output:
```
row.names(rna_dna_path_unstratified)[1:10]
```
```
  [1] "1CMET2-PWY: N10-formyl-tetrahydrofolate biosynthesis"                                                          
  [2] "3-HYDROXYPHENYLACETATE-DEGRADATION-PWY: 4-hydroxyphenylacetate degradation"                                    
  [3] "7ALPHADEHYDROX-PWY: cholate degradation (bacteria, anaerobic)"                                                 
  [4] "AEROBACTINSYN-PWY: aerobactin biosynthesis"                                                                    
  [5] "ALL-CHORISMATE-PWY: superpathway of chorismate metabolism"                                                     
  [6] "ALLANTOINDEG-PWY: superpathway of allantoin degradation in yeast"                                              
  [7] "ANAEROFRUCAT-PWY: homolactic fermentation"                                                                     
  [8] "ANAGLYCOLYSIS-PWY: glycolysis III (from glucose)"                                                              
  [9] "ARG+POLYAMINE-SYN: superpathway of arginine and polyamine biosynthesis"                                        
 [10] "ARGDEG-PWY: superpathway of L-arginine, putrescine, and 4-aminobutanoate degradation"    
```

Remove temp files to clear up space:
```
rm(tmp.ind)
```

Now that we have sucessfully unstratified `rna_dna_path`, let's check the sample sums (colSums):
```
colSums(rna_dna_path_unstratified)
```
```
   MSM5LLDI    MSM9VZL5    PSM6XBRK    PSMA265X    MSM79H8D    PSM7J4EF 
  208.80878   306.06509   175.91123   403.28804   297.47223   184.13486 
...
```


* Question: Why aren't these samples summing to ~1 like all of the datasets before?


Filter for beta diversity:

Only keep RNA/DNA pathways that passed filtering for DNA pathways:

Check the dimensions of RNA/DNA pathways pre-filtering:
```
dim(rna_dna_path_unstratified)
```
```
[1] 472  12
```
Filter:
```
rna_dna_path_unstratified_filt = subset(rna_dna_path_unstratified, row.names(rna_dna_path_unstratified) %in% names(dna_path_unstratified_filt))
```
Check the dimensions of RNA/DNA pathways post-filtering:
```
dim(rna_dna_path_unstratified_filt)
```
```
[1] 309  12
```
Check the dimensions of dna pathways post-filtering to make sure they line up:
```
dim(dna_path_unstratified_filt)
```
```
[1]  96 309
```
Yes, they have the same number of pathways.


Let's transpose the dataframes for easier use downstream, making the rows be the samples just like in the metadata:
```
rna_dna_path_unstratified_filt = data.frame(t(rna_dna_path_unstratified_filt), check.names = F)
rna_dna_path_unstratified = data.frame(t(rna_dna_path_unstratified), check.names = F)
```

Check to make sure the sample names are in the same order in both the metadata and the RNA/DNA pathways dataframe:
```
row.names(metadata_rna_dna) == row.names(rna_dna_path_unstratified_filt)
```
```
 [1] FALSE FALSE FALSE FALSE  TRUE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE
 ```
They don't match, so let's match them up.

 ```
rna_dna_path_unstratified_filt = rna_dna_path_unstratified_filt[match(row.names(metadata_rna_dna), row.names(rna_dna_path_unstratified_filt)),]
```
```
row.names(metadata_rna_dna) == row.names(rna_dna_path_unstratified_filt)
```
```
 [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
```
Good to go now.

>NOTE: These formatted files are also located in the Tutorials/highdimtesting directory of the bioBakery image. To work with them from there just assign them in R with the following code:
>```
>metadata_rna_dna = read.csv(file = "metadata_rna_dna.csv", header = T, row.names = 1, check.names = FALSE)
>rna_dna_path_unstratified_filt = read.csv(file = "rna_dna_path_unstratified_filt.csv", header = T, row.names = 1, check.names = FALSE)
>```

log transform the RNA/DNA ratio:
```
rna_dna_path_unstratified_filt_log = log2(rna_dna_path_unstratified_filt + 1)
```


### Beta Diversity PERMANOVA tests using Euclidean distances

Calculate Euclidean distances:
```
euclidean_rna_dna_path_unstratified = vegdist(rna_dna_path_unstratified_filt_log, method = "euclidean")
```

* Question: Why are we using Euclidean distances on this dataset as opposed to Bray-Curtis dissimilarities that we used on all of the other feature tables?


#### Univariable

Can do this with a for loop:
```
for (col in names(metadata_rna_dna)){
  adonis_univ <- adonis(as.formula(paste("euclidean_rna_dna_path_unstratified ~ ", col)), data = metadata_rna_dna)
  print(adonis_univ)
}
```
```
Call:
adonis(formula = as.formula(paste("euclidean_rna_dna_path_unstratified ~ ",      col)), data = metadata_rna_dna) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
site_name  2    365.27  182.63  1.5789 0.25974  0.022 *
Residuals  9   1041.03  115.67         0.74026         
Total     11   1406.29                 1.00000         
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Call:
adonis(formula = as.formula(paste("euclidean_rna_dna_path_unstratified ~ ",      col)), data = metadata_rna_dna) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
sex        1    148.92  148.91  1.1843 0.10589  0.204
Residuals 10   1257.38  125.74         0.89411       
Total     11   1406.29                 1.00000       
Error in `contrasts<-`(`*tmp*`, value = contr.funs[1 + isOF[nn]]) : 
  contrasts can be applied only to factors with 2 or more levels
```

* Question: We got an error, what could be the cause?
    * Hint: Let's investigate the metadata for this set of samples again.


```
metadata_rna_dna
```
```
                 site_name    sex  race consent_age diagnosis
HSM6XRQB        Cincinnati Female White          16        CD
HSM7J4LP        Cincinnati Female White          12        CD
MSM5LLDI               MGH Female White          30        CD
MSM6J2Q1               MGH   Male White          28    nonIBD
MSM79H8D               MGH   Male White          61    nonIBD
MSM9VZL5               MGH   Male White          21    nonIBD
MSM9VZMM               MGH Female White          23    nonIBD
PSM6XBRK    MGH Pediatrics   Male White          16        CD
PSM6XBRK_TR MGH Pediatrics   Male White          16        CD
PSM6XBSE    MGH Pediatrics Female White          16        UC
PSM7J4EF    MGH Pediatrics   Male White          15        CD
PSMA265X    MGH Pediatrics Female White          16        UC
```

We only have one race now in this dataset, so let's remove that column from the metadata in the formula call. This wouldn't be a problem in the one-by-one method below, but is in the loop that goes through each column one at a time. Remember how we did this before for the taxonomy data in the alpha diversity models?
```
for (col in names(metadata_rna_dna)[-3]){
  adonis_univ <- adonis(as.formula(paste("euclidean_rna_dna_path_unstratified ~ ", col)), data = metadata_rna_dna[-3])
  print(adonis_univ)
}
```
```
Call:
adonis(formula = as.formula(paste("euclidean_rna_dna_path_unstratified ~ ",      col)), data = metadata_rna_dna[-3]) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
site_name  2    365.27  182.63  1.5789 0.25974  0.017 *
Residuals  9   1041.03  115.67         0.74026         
Total     11   1406.29                 1.00000         
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Call:
adonis(formula = as.formula(paste("euclidean_rna_dna_path_unstratified ~ ",      col)), data = metadata_rna_dna[-3]) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
sex        1    148.92  148.91  1.1843 0.10589  0.217
Residuals 10   1257.38  125.74         0.89411       
Total     11   1406.29                 1.00000       

Call:
adonis(formula = as.formula(paste("euclidean_rna_dna_path_unstratified ~ ",      col)), data = metadata_rna_dna[-3]) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
consent_age  1    152.54  152.54  1.2167 0.10847  0.186
Residuals   10   1253.75  125.38         0.89153       
Total       11   1406.29                 1.00000       

Call:
adonis(formula = as.formula(paste("euclidean_rna_dna_path_unstratified ~ ",      col)), data = metadata_rna_dna[-3]) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
diagnosis  2    478.92  239.46  2.3239 0.34056  0.001 ***
Residuals  9    927.37  103.04         0.65944           
Total     11   1406.29                 1.00000           
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

Alternatively, can do it one-by-one like above in the taxonomy section.


#### Multivariable

Can do it without being verbose:
```
adonis_multi_rna_dna_path = adonis(euclidean_rna_dna_path_unstratified ~ ., data = metadata_rna_dna[-3])
adonis_multi_rna_dna_path
```
```
Call:
adonis(formula = euclidean_rna_dna_path_unstratified ~ ., data = metadata_rna_dna[-3]) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
site_name    2    365.27 182.633 1.78876 0.25974  0.011 *
sex          1    139.76 139.756 1.36880 0.09938  0.106  
consent_age  1     87.34  87.337 0.85541 0.06210  0.615  
diagnosis    2    303.43 151.715 1.48594 0.21577  0.050 *
Residuals    5    510.50 102.100         0.36301         
Total       11   1406.29                 1.00000         
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

Diagnosis goes from being very significant in the univariable model to nearly non-significant in the multivariable.

Alternatively, can write out each variable in the model like above in the taxonomy section.


#### Pairwise Euclidean comparisons
Pairwise comparisons of diagnosis in a univariable model

```
euclidean_rna_dna_path_unstratified_mat = as.matrix(euclidean_rna_dna_path_unstratified)
pairwise.adonis(euclidean_rna_dna_path_unstratified_mat, factors = metadata_rna_dna$diagnosis, p.adjust.m = "BH")
```
```
'nperm' >= set of all permutations: complete enumeration.
Set of permutations < 'minperm'. Generating entire set.
         pairs Df  SumsOfSqs  F.Model        R2   p.value p.adjusted sig
1 CD vs nonIBD  1 0.03970044 5.395477 0.4027835 0.0030000  0.0090000   *
2     CD vs UC  1 0.01941921 2.108531 0.2600386 0.0660000  0.0990000    
3 nonIBD vs UC  1 0.01923794 1.947107 0.3274041 0.1333333  0.1333333    
```


Now let's look at all of the PERMANOVA results together:
```
adonis_multi_tax
```
```
Call:
adonis(formula = bray_species ~ ., data = metadata) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
site_name    4    1.6651 0.41628  1.4195 0.05767  0.015 *
sex          1    0.4848 0.48479  1.6531 0.01679  0.038 *
race         3    0.8897 0.29658  1.0113 0.03081  0.433  
consent_age  1    0.3545 0.35448  1.2088 0.01228  0.240  
diagnosis    2    0.8470 0.42349  1.4441 0.02933  0.046 *
Residuals   84   24.6335 0.29326         0.85312         
Total       95   28.8747                 1.00000         
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
```
adonis_multi_dna_path
```
```
Call:
adonis(formula = bray_dna_path_unstratified ~ ., data = metadata) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
site_name    4    0.1763 0.044076  1.1105 0.04646  0.273
sex          1    0.0311 0.031101  0.7836 0.00820  0.616
race         3    0.1243 0.041425  1.0437 0.03275  0.394
consent_age  1    0.0293 0.029283  0.7378 0.00772  0.620
diagnosis    2    0.1001 0.050034  1.2606 0.02637  0.222
Residuals   84    3.3339 0.039690         0.87852       
Total       95    3.7950                  1.00000       
```
```
adonis_multi_rna_path
```
```
Call:
adonis(formula = bray_rna_path_unstratified ~ ., data = metadata_rna) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
site_name    4   0.33933 0.084834 0.94530 0.15419  0.489
sex          1   0.05461 0.054609 0.60851 0.02481  0.888
race         3   0.16446 0.054820 0.61086 0.07473  0.819
consent_age  1   0.06739 0.067385 0.75088 0.03062  0.558
diagnosis    2   0.13903 0.069515 0.77461 0.06318  0.811
Residuals   16   1.43587 0.089742         0.65246       
Total       27   2.20069                  1.00000       
```
```
adonis_multi_rna_dna_path
```
```
Call:
adonis(formula = euclidean_rna_dna_path_unstratified ~ ., data = metadata_rna_dna[-3]) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
site_name    2    365.27 182.633 1.78876 0.25974  0.011 *
sex          1    139.76 139.756 1.36880 0.09938  0.106  
consent_age  1     87.34  87.337 0.85541 0.06210  0.615  
diagnosis    2    303.43 151.715 1.48594 0.21577  0.050 *
Residuals    5    510.50 102.100         0.36301         
Total       11   1406.29                 1.00000         
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

* Question: What are some possible reasons the taxonomy data gives the only significant results (aside from the RNA/DNA ratio data that could be misleading due to low n)?


# Mantel tests

To be continued...
