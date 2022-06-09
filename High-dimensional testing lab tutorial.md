# High-dimensional testing tutorial

- [Omnibus tests on data from the HMP2 baseline only data](#omnibus-tests-on-data-from-the-hmp2-baseline-only-data)
  * [MGX taxonomy](#mgx-taxonomy)
    + [Feature table and metadata table creation and formatting](#feature-table-and-metadata-table-creation-and-formatting)
    + [Alpha diversity](#alpha-diversity)
      - [Univariable](#univariable)
      - [Multivariable](#multivariable)
    + [Beta Diversity](#beta-diversity)
      - [Univariable](#univariable-1)
      - [Multivariable](#multivariable-1)

# Omnibus tests on data from the HMP2 baseline only data

In this tutorial, we are going to take data from the second phase of the Human Microbiome Project, and look at some common exploratory analyses which search for evidence of general microbiome associations with cohort metadata. These tests are frequently the first we perform on a new dataset (after quality control) and can give a good idea of things such as the overall strength of an intervention or the degree to which covariates might confound the main microbiome outcomes of interest.

We are starting with the baseline timepoint, metagenome-derived species data, since species tables (+ associated metadata) are what you will likely most commonly encounter, and longitudinal designs require more advanced statistical techniques. 

Start R and set the working directory. In the bioBakery VM, as follows:
```
R
setwd("~/Tutorials/highdimtesting")
```

Load the R packages needed:
```
library(vegan)
```

---
* **What are omnibus tests and how do they differ from featurewise tests?**
* **How might these complement each other?**
<details>
 <summary><b>Answer:</b></summary>
 
 Omnibus tests are global tests, here of the microbial community diversity. Especially in smaller studies, it's possible that no individual features reach significance, but some microbiome association can nonetheless be statistically shown.
</details>

---

## MGX taxonomy

### Feature table and metadata table creation and formatting

Download the MGX taxonomy relative abundance data and put it into the data directory:
```
download.file("https://raw.githubusercontent.com/biobakery/omnibus-and-maaslin2-rscripts-and-hmp2-data/master/taxonomic_profiles_pcl_week0.csv", "./Data/taxonomic_profiles_pcl_week0.csv")
```

Read the taxonomic relative abundance data into your R environment:
```
tax = read.csv(file = "./Data/taxonomic_profiles_pcl_week0.csv", header = T, row.names = 1, check.names = F)
```

Take a look at `tax`:
```
tax[1:10,1:10]
```

---
* **Why not use head(tax)?** Try it if you are unsure.

<details>
 <summary><b>Answer:</b></summary>
 
 `head` cuts off at a certain number of rows, but not columns.
</details>

---

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

---
* **What are the advantages and disadvantages of storing and accessing data like this?** 

<details>
 <summary><b>Answer:</b></summary>
 
 A single file can be easier to deal with and is required for some tools, but accessing data with static indicies can cause problems if the file is changed (e.g. if another metadata column is added then tax[1:5] is not correct). 
</details>

---

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

---
* **What is the main drawback of keeping NA values for discrete variables?** 
* **Is there a case where this is totally justified?**

<details>
 <summary><b>Answer:</b></summary>
 
 Adding another category means increasing the degrees of freedom in downstream statistical tests. You may however want to specifically investigate if missing values are meaningful (e.g. failure to report a value can indicate study attrition which may be caused by an adverse reaction).
</details>

---

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

Verify that strains were removed and skim the output to make sure nothing unexpected happened:
```
grep("\\|t__", rownames(tmp), invert = F)
row.names(tmp)[1:10]
```
```
integer(0)
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

---
* **What will happen if there is no strain level information in the starting file?**

<details>
 <summary><b>Answer:</b></summary>
 
 There's no issue, but it's always good to check that your code can handle different situations. MetaPhlAn doesn't alway return strain stratifications, so it's is particularly important in this case.
</details>

---

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
colSums(species)[1:6]
```
```
   CSM67UH7    CSM79HHW    HSM67VDT    HSM6XRQB    HSM6XRR3    HSM7J4LP 
  0.9999752   1.0000004   1.0000005   0.9993355   0.9996256   1.0000003 
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

---
* **What exactly is this command doing? What does 0.0001 represent? 0.1?** Note that you can break the command apart to poke at the output, e.g. run `apply(species, 1, function(x) sum(x > 0.0001)` first.

<details>
 <summary><b>Answer:</b></summary>
 
 This filters out species with less than 0.0001 relative abundance in 10% of samples.
</details>

---

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

Make sure the sample names are in the same order in both the metadata and the species dataframes:
```
all(row.names(metadata) == row.names(species_filt))
all(row.names(metadata) == row.names(species))
```

---
* **Why perform this type of filtering? Are there situations were it might not be appropriate?**

<details>
 <summary><b>Answer:</b></summary>
 
 Low prevalence/abundance features are often not meaningful and can reduce power of some tests. However, one needs to be certain this isn't removing real biological variation. Alpha diversity is sensitive to filtering (so we generally don't) and some datasets have important low prevalence taxa or outlier samples.
</details>

---

>NOTE: These formatted files are also located in the Tutorials/highdimtesting directory of the bioBakery image. To work with them from there assign them in R with the following code:
>```
>metadata = read.csv(file = "metadata.csv", header = T, row.names = 1, check.names = F)
>species = read.csv(file = "species.csv", header = T, row.names = 1, check.names = F)
>species_filt = read.csv(file = "species_filt.csv", header = T, row.names = 1, check.names = F)
>```

### Alpha diversity

---
* **What is alpha diversity, and why is it important?**

<details>
 <summary><b>Answer:</b></summary>
 
 Alpha diversity summarizes the structure of an ecological community with respect to its richness (roughly, number of features), evenness (roughly, distribution of features), or both. Many perturbations to the microbiome affect this.
</details>

---

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
         diversity(species, index = "shannon")
CSM67UH7                             1.7906860
CSM79HHW                             2.5342475
HSM67VDT                             2.2925910
HSM6XRQB                             2.5677991
HSM6XRR3                             2.2772366
HSM7J4LP                             0.3959699
```

Rename the columns to shorter names:
```
names(alpha_div) = c("Shannon")
```

Collate with the metadata, so that a single dataframe can be provided to upcoming statistical tests:
```
alpha_meta <- merge(alpha_div, metadata, by = "row.names")
```
Note that the `merge` function makes the old row names into new column. 
This might not be what you expect, but you can see that this is intended with `?merge`.

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
            Shannon    site_name    sex  race consent_age diagnosis
CSM5FZ3N_P 1.741766 Cedars-Sinai Female White          43        CD
CSM5FZ3T_P 1.634631 Cedars-Sinai Female White          76        CD
CSM5FZ4A_P 2.117305 Cedars-Sinai Female White          47        UC
CSM5MCTZ_P 2.665711 Cedars-Sinai   Male White          32        UC
CSM5MCU4_P 1.119533 Cedars-Sinai Female White          53        CD
CSM5MCV1_P 1.478073 Cedars-Sinai Female White          20        CD
```

#### Univariable

Alpha diversity is relatively well-behaved and generates a single value per sample, allowing for a number of statistical tests.
Linear models are straightforward and can accommodate a large variety of experimental designs, so we will go with that.

Run a linear model on diagnosis and return as an ANOVA table:
```
anova(lm(Shannon ~ diagnosis, data = alpha_meta))
```
```
Analysis of Variance Table

Response: Shannon
          Df  Sum Sq Mean Sq F value Pr(>F)
diagnosis  2  0.4046 0.20231  0.6692 0.5146
Residuals 93 28.1163 0.30233               
```

Can do this with a for loop to quickly get results for each metadata variable:
```
for (col in names(metadata)) {
  alpha_shannon_univ = anova(lm(as.formula(paste("Shannon ~ ", col)), data = alpha_meta)) 
  print(col)
  print(alpha_shannon_univ)
}
```
```
[1] "site_name"
Analysis of Variance Table

Response: Shannon
          Df  Sum Sq Mean Sq F value Pr(>F)
site_name  4  1.6732 0.41831  1.4179 0.2344
Residuals 91 26.8477 0.29503               
[1] "sex"
Analysis of Variance Table

Response: Shannon
          Df  Sum Sq Mean Sq F value  Pr(>F)  
sex        1  0.9209 0.92088  3.1363 0.07981 .
Residuals 94 27.6001 0.29362                  
---
Signif. codes:  0 `***` 0.001 `**` 0.01 `*` 0.05 `.` 0.1 ` ` 1
[1] "race"
Analysis of Variance Table

Response: Shannon
          Df  Sum Sq Mean Sq F value Pr(>F)
race       3  0.9266 0.30886  1.0297 0.3832
Residuals 92 27.5944 0.29994               
[1] "consent_age"
Analysis of Variance Table

Response: Shannon
            Df  Sum Sq  Mean Sq F value Pr(>F)
consent_age  1  0.0105 0.010538  0.0347 0.8525
Residuals   94 28.5104 0.303302               
[1] "diagnosis"
Analysis of Variance Table

Response: Shannon
          Df  Sum Sq Mean Sq F value Pr(>F)
diagnosis  2  0.4046 0.20231  0.6692 0.5146
Residuals 93 28.1163 0.30233               
```

#### Multivariable

---
* **What are univariable and multivariable tests, and what benefits do both provide?**

<details>
 <summary><b>Answer:</b></summary>
 
 Univariable tests have only one independent variable and multivariable have more. Note that multivariate formally refers to multiple dependent variables, though these terms are often used interchangeably. 
</details>

---

Try a model with the age and diagnosis:
```
anova(lm(Shannon ~ diagnosis + sex, data = alpha_meta))
```
```
Analysis of Variance Table

Response: Shannon
          Df  Sum Sq Mean Sq F value Pr(>F)
diagnosis  2  0.4046 0.20231  0.6808 0.5087
sex        1  0.7781 0.77807  2.6184 0.1091
Residuals 92 27.3383 0.29716  
```

If all variables in the file are used, there's a nice shorthand in R:
```
anova(lm(Shannon ~ ., data = alpha_meta))
```
```
Analysis of Variance Table

Response: Shannon
            Df  Sum Sq Mean Sq F value Pr(>F)
site_name    4  1.6732 0.41831  1.4118 0.2371
sex          1  0.4001 0.40012  1.3504 0.2485
race         3  0.6939 0.23130  0.7807 0.5080
consent_age  1  0.7356 0.73559  2.4826 0.1189
diagnosis    2  0.1294 0.06472  0.2184 0.8042
Residuals   84 24.8887 0.29629               
```

---
* **This can be expanded to include random effects with the `lme4` package. Is there anything in this demo data that might indicate a mixed model would be useful?**

<details>
 <summary><b>Answer:</b></summary>
 
 `site_name` is the location where the sample was collected and is likely better considered a random effect, since we aren't necessarily interested in how the sites vary but rather controlling for the overall effect of sampling site (of theoretically many).
</details>

---

### Beta Diversity

---
* **What is beta diversity, and why is it important?**

<details>
 <summary><b>Answer:</b></summary>
 
 Beta diversity quantifies how similar (or dissimilar) pairs of samples are in composition. This can be used to create a map of distances between samples (e.g. PCoA plot) and the significance of grouping or gradients of samples can be tested.
</details>

---

Using the filtered species table, calculate Bray-Curtis dissimilarity:
```
bray_species = vegdist(species_filt, method = "bray")
```
Bray-Curtis is commonly used for microbiome data because it incorporates presence vs. absence and relative abundance of features when determining dissimilarities (i.e. is weighted). Another common index for species data is Weighted UniFrac, which additionally includes phylogenetic information. 

#### Univariable

Unlike alpha diversity, beta diversity is not defined on a single sample level, and therefore, many common statistical tests are not usable. 

---
* **Why bother with beta diversity then?**

<details>
 <summary><b>Answer:</b></summary>
 
 While one could take a difference of alpha diversity between samples and get a number, this is not particularly robust or meaningful. Alpha diversity trends are also in general more sensitive to difficult to control or arbitrary factors, such as sampling depth and technical noise.
</details>

---

PERMANOVA is the most common omnibus beta diversity test and is implemented as `adonis2` by the vegan package.
Note that `adonis` is now officially deprecated and older code that uses this function may act strangely. 

Run a PERMANOVA and print the results:
```
set.seed(123)
adonis2(bray_species ~ diagnosis, data = metadata)
```
```
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = bray_species ~ diagnosis, data = metadata)
          Df SumOfSqs      R2      F Pr(>F)  
diagnosis  2   0.9399 0.03255 1.5646   0.02 *
Residual  93  27.9347 0.96745                
Total     95  28.8747 1.00000                
---
Signif. codes:  `***` 0.001 `**` 0.01 `*` 0.05 `.` 0.1 ` ` 1
```
As you can see, the results format looks quite similar to the more familiar ANOVA test.

---
* **Why is a seed set before calling `adonis2`?** Try running it again without resetting the seed.

<details>
 <summary><b>Answer:</b></summary>
 
 PERMANOVA determines p-values via permutation tests, so unless a seed is set, the values will be different each time you run it.
</details>

---

Can do this with a for loop as well:
```
for (col in names(metadata)) {
  set.seed(123)
  adonis_univ <- adonis2(as.formula(paste("bray_species ~ ", col)), data = metadata)
  print(adonis_univ)
}
```
```
[1] "site_name"
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = as.formula(paste("bray_species ~ ", col)), data = metadata)
          Df SumOfSqs      R2      F Pr(>F)  
site_name  4   1.6651 0.05767 1.3922  0.021 *
Residual  91  27.2095 0.94233                
Total     95  28.8747 1.00000                
---
Signif. codes:  0 `***` 0.001 `**` 0.01 `*` 0.05 `.` 0.1 ` ` 1
[1] "sex"
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = as.formula(paste("bray_species ~ ", col)), data = metadata)
         Df SumOfSqs      R2      F Pr(>F)  
sex       1   0.5925 0.02052 1.9694  0.013 *
Residual 94  28.2821 0.97948                
Total    95  28.8747 1.00000                
---
Signif. codes:  0 `***` 0.001 `**` 0.01 `*` 0.05 `.` 0.1 ` ` 1
[1] "race"
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = as.formula(paste("bray_species ~ ", col)), data = metadata)
         Df SumOfSqs      R2      F Pr(>F)
race      3   0.7778 0.02694 0.8489  0.771
Residual 92  28.0969 0.97306              
Total    95  28.8747 1.00000              
[1] "consent_age"
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = as.formula(paste("bray_species ~ ", col)), data = metadata)
            Df SumOfSqs      R2      F Pr(>F)  
consent_age  1   0.4683 0.01622 1.5497  0.067 .
Residual    94  28.4063 0.98378                
Total       95  28.8747 1.00000                
---
Signif. codes:  0 `***` 0.001 `**` 0.01 `*` 0.05 `.` 0.1 ` ` 1
[1] "diagnosis"
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = as.formula(paste("bray_species ~ ", col)), data = metadata)
          Df SumOfSqs      R2      F Pr(>F)  
diagnosis  2   0.9399 0.03255 1.5646   0.02 *
Residual  93  27.9347 0.96745                
Total     95  28.8747 1.00000                
---
Signif. codes:  0 `***` 0.001 `**` 0.01 `*` 0.05 `.` 0.1 ` ` 1
```

---
* **What does the R2 value tell us?** 
* **Biologically, when might you expect a large R2? How about a small but significant R2?**

<details>
 <summary><b>Answer:</b></summary>
 
 R2 is proportion of the variance for a dependent variable that's explained by an independent variable. You might expect a large R2 for a distinct perturbation like recent antibiotic treatment and a small but significant R2 for a complex population level trend, such as urbanicity in childhood. 
</details>

---

#### Multivariable

Try a model with age and diagnosis:
```
set.seed(123)
adonis2(bray_species ~ consent_age + diagnosis, data = metadata)
```
```
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = bray_species ~ consent_age + diagnosis, data = metadata)
            Df SumOfSqs      R2      F Pr(>F)  
consent_age  1   0.4683 0.01622 1.5696  0.060 .
diagnosis    2   0.9569 0.03314 1.6035  0.015 *
Residual    92  27.4495 0.95064                
Total       95  28.8747 1.00000                
---
Signif. codes:  0 `***` 0.001 `**` 0.01 `*` 0.05 `.` 0.1 ` ` 1
```

PERMANOVA is by default sequential, but it's possible to do a marginal test:
```
set.seed(123)
adonis2(bray_species ~ consent_age + diagnosis, data = metadata, by = "margin")
```
```
Permutation test for adonis under reduced model
Marginal effects of terms
Permutation: free
Number of permutations: 999

adonis2(formula = bray_species ~ consent_age + diagnosis, data = metadata, by = "margin")
            Df SumOfSqs      R2      F Pr(>F)  
consent_age  1   0.4853 0.01681 1.6264  0.048 *
diagnosis    2   0.9569 0.03314 1.6035  0.015 *
Residual    92  27.4495 0.95064                
Total       95  28.8747 1.00000                
---
Signif. codes:  0 `***` 0.001 `**` 0.01 `*` 0.05 `.` 0.1 ` ` 1
```
Note that this can get very computationally expensive as more terms are added.

---
* **What is the difference between these two models?** 

<details>
 <summary><b>Answer:</b></summary>
 
 The marginal model is agnostic to the ordering of terms. For PERMANOVA, it can be generated by hand by running the sequential model twice (see that the second term, diagnosis, has the same output both times).  The trade-off is that the R2 do not necessarily add up to 1.  

</details>

---

Can use the same shorthand as with linear models to include all variables:
```
set.seed(123)
adonis2(bray_species ~ ., data = metadata)
```
```
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = bray_species ~ ., data = metadata)
            Df SumOfSqs      R2      F Pr(>F)  
site_name    4   1.6651 0.05767 1.4195  0.018 *
sex          1   0.4848 0.01679 1.6531  0.044 *
race         3   0.8897 0.03081 1.0113  0.438  
consent_age  1   0.3545 0.01228 1.2088  0.230  
diagnosis    2   0.8470 0.02933 1.4441  0.051 .
Residual    84  24.6335 0.85312                
Total       95  28.8747 1.00000                
---
Signif. codes:  0 `***` 0.001 `**` 0.01 `*` 0.05 `.` 0.1 ` ` 1
```

Interaction terms are also possible:
```
set.seed(123)
adonis2(bray_species ~ diagnosis * sex, data = metadata)
```
```
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = bray_species ~ diagnosis * sex, data = metadata)
              Df SumOfSqs      R2      F Pr(>F)  
diagnosis      2   0.9399 0.03255 1.5909  0.015 *
sex            1   0.5614 0.01944 1.9004  0.013 *
diagnosis:sex  2   0.7874 0.02727 1.3327  0.059 .
Residual      90  26.5860 0.92074                
Total         95  28.8747 1.00000                
---
Signif. codes:  0 `***` 0.001 `**` 0.01 `*` 0.05 `.` 0.1 ` ` 1
```

---
* **Is this result reasonable?** 
* **What factors might complicate this relationship?** 

<details>
 <summary><b>Answer:</b></summary>
 
 There are known epidemiological incidence and prevalence sex differences in IBD, but they are not straightforward and are confounded by other factors such as age and type of IBD (UC vs CD). Also, some caution should be taken when using and interpreting interaction effects in PERMANOVA models since these are not as developed or rigorously tested as in ANOVA.
</details>

---
