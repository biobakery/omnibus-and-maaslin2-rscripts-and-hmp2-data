# High-dimensional testing lab tutorial

# Omnibus tests on mgx and mtx data from the HMP2. Baseline only data.

Create new directories and set the working directory:
```
dir.create("R_Metagenomic_Statistics")
dir.create("R_Metagenomic_Statistics/Data")
setwd("R_Metagenomic_Statistics")
```

Load the R packages needed:
```
library(vegan)
library(plyr)
```

* Question: What are omnibus tests and how do they differ from featurewise tests? What types of plots do these omnibus tests complement?

## MGX taxonomy

### Feature table and metadata table creation and formatting

Download the mgx taxonomy relative abundance data and put it into the data directory:
```
download.file("https://raw.githubusercontent.com/biobakery/omnibus-and-maaslin2-rscripts-and-hmp2-data/master/taxonomic_profiles_pcl_week0.csv", "./Data/taxonomic_profiles_pcl_week0.csv")
```
Read the taxonomic relative abundance data into your R environment:
```
tax = read.csv(file = "./Data/taxonomic_profiles_pcl_week0.csv", header = T, row.names = 1, check.names = FALSE)
```

Check out `tax`:

Check the dimensions:
```
dim(tax)

[1]   96 1484
```

Extract the metadata from `tax`:
```
metadata = data.frame(tax[1:5])
```

Check out `metadata`:

Check the dimensions:
```
dim (metadata)

[1] 96  5
```
Check the output:
```
head(metadata)

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

'data.frame':	96 obs. of  5 variables:
 $ site_name  : Factor w/ 5 levels "Cedars-Sinai",..: 1 1 2 2 2 2 2 2 4 4 ...
 $ sex        : Factor w/ 2 levels "Female","Male": 2 1 2 1 2 1 2 1 1 2 ...
 $ race       : Factor w/ 4 levels "Black or African American",..: 4 4 4 4 4 4 4 4 4 4 ...
 $ consent_age: int  69 36 11 16 13 12 16 14 30 57 ...
 $ diagnosis  : Factor w/ 3 levels "CD","nonIBD",..: 2 1 2 1 1 1 3 2 1 2 ...
 ```

 Check for NAs in `metadata` that will later cause issues with the PERMANOVAs:
 ```
count(is.na(metadata))

  x.site_name x.sex x.race x.consent_age x.diagnosis freq
1       FALSE FALSE  FALSE         FALSE       FALSE   90
2       FALSE FALSE  FALSE          TRUE       FALSE    6
```
Age has 6 NA.

If this was a discrete variable we could just classify the NAs as Unknown and keep them in the model, but since age is a continuous variable typically we would either remove those from the data or impute the median. 

In this case let's impute the median in order to keep samples:

Check out unique ages:
```
unique(metadata$consent_age)

[1] 69 36 11 16 13 12 14 30 57 28 61 29 51 21 25 23 17 15 43 47 NA  7 10  8 37 44 26 32  9 40 74 19 55 50 45  6 56 38 76 24 41 53
```
Notice the NA within the numbers.

```
metadata$consent_age[is.na(metadata$consent_age)] = median(metadata$consent_age, na.rm = T)
```
Check the output:
```
unique(metadata$consent_age)

[1] 69 36 11 16 13 12 14 30 57 28 61 29 51 21 25 23 17 15 43 47 20  7 10  8 37 44 26 32  9 40 74 19 55 50 45  6 56 38 76 24 41 53
```
NA has been replaced by 20; good to go.


Extract species data from `tax` and transpose the dataframe:
```
species = data.frame(t(tax[6:ncol(tax)]))
```

Check out `species`:

Check out the structure:
```
str(species)

'data.frame':	1479 obs. of  96 variables:
 $ CSM67UH7   : num  0 0 0 0 0 ...
 $ CSM79HHW   : num  0 0 0 0 0 0 0 0 0 1 ...
 $ HSM67VDT   : num  0 0 0 0 0 0 0 0 0 1 ...
 $ HSM6XRQB   : num  0 0 0 0 0 0 0 0 0 1 ...
 $ HSM6XRR3   : num  0 0 0 0 0 0 0 0 0 1 ...
 $ HSM7J4LP   : num  0 0 0 0 0 0 0 0 0 1 ...
 ...
```
Everything is numeric; good to go.

Check the output:
```
species[1:8,1:2]

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
As we can see here, `species` currently has all taxonomic levels for each microbe, so now we need to only keep the species level information.

grep the rows that do not include strain stratifications and store them in a new temporary vector:
```
tmp.ind = grep("\\|t__", rownames(species), invert = T)
```

Create a new dataframe with only those row numbers in `tmp.ind`:
```
tmp = species[tmp.ind,]
```
Check the output:
```
row.names(tmp)

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
 ...
```

Now that we have removed the strain stratifications, let's select only the species level stratifications, removing all taxonomic levels before it:

```
tmp.ind = grep("\\|s__", rownames(tmp))
```

Create a new dataframe with only those row numbers in `tmp.ind`:
```
species = tmp[tmp.ind,]
```

Check the output to make sure that we only have species level stratifications:
```
row.names(species)

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
 ...
```

Remove temp files to clear up space:
```
rm(tmp,tmp.ind)
```

Now we have only species level information in `species`, but let's make the names shorter by trimming off all other taxonomic information that predates species:
```
rownames(species) = gsub(".*\\|", "", rownames(species))
```
Check the output:
```
row.names(species)

  [1] "s__Methanobrevibacter_smithii"                 "s__Methanobrevibacter_unclassified"            "s__Granulicella_unclassified"                 
  [4] "s__Terriglobus_unclassified"                   "s__Actinobaculum_schaalii"                     "s__Actinobaculum_unclassified"                
```

Now that we have sucessfully truncated the names in `species`, let's check the sample sums (colSums) to make sure they are in proportion format (0-1) and are all ~1:
```
colSums(species)

   CSM67UH7    CSM79HHW    HSM67VDT    HSM6XRQB    HSM6XRR3    HSM7J4LP 
  0.9999752   1.0000004   1.0000005   0.9993355   0.9996256   1.0000003 
...
```

Filter for beta diversity (we will keep species as is for alpha diversity):

Check the dimensions of species pre-filtering:
```
dim(species)

[1] 572  96
```
Filter:
```
species_filt = species[apply(species, 1, function(x) sum(x > 0.0001) > 0.1 * ncol(species)), ]
```
Check the dimensions of species post-filtering:
```
dim(species_filt)

[1] 115  96
```

Let's transpose the dataframes for easier use downstream, making the rows be the samples just like in the metadata:
```
species_filt = data.frame(t(species_filt), check.names = F)
species = data.frame(t(species), check.names = F)
```


### Alpha diversity

Using unfiltered species data, we will calculate a few different alpha diveristy metrics:
```
shannon = data.frame(diversity(species, index = "shannon"))
simpson = data.frame(diversity(species, index = "simpson"))
invsimpson = data.frame(diversity(species, index = "invsimpson"))
```

* Question: What is alpha diversity and why is it important? 

Merge all into one dataframe:
```
alphadiv = cbind(shannon, simpson, invsimpson)
```
Check the output:
```
head(alphadiv)

         diversity.species..index....shannon.. diversity.species..index....simpson.. diversity.species..index....invsimpson..
CSM67UH7                             1.7906860                             0.6360532                                 2.747654
CSM79HHW                             2.5342475                             0.8963525                                 9.648086
HSM67VDT                             2.2925910                             0.7982095                                 4.955635
HSM6XRQB                             2.5677991                             0.8790000                                 8.264464
HSM6XRR3                             2.2772366                             0.8455911                                 6.476311
HSM7J4LP                             0.3959699                             0.1864280                                 1.229147
```

Rename the columns to shorter names:
```
names(alphadiv) = c("Shannon", "Simpson", "InvSimpson")
```
Check the output:
```
head(alphadiv)

           Shannon   Simpson InvSimpson
CSM67UH7 1.7906860 0.6360532   2.747654
CSM79HHW 2.5342475 0.8963525   9.648086
HSM67VDT 2.2925910 0.7982095   4.955635
HSM6XRQB 2.5677991 0.8790000   8.264464
HSM6XRR3 2.2772366 0.8455911   6.476311
HSM7J4LP 0.3959699 0.1864280   1.229147
```

Append metadata:
```
alpha_met_df <- merge(alphadiv, metadata, by = "row.names")
```

Make the samples the row names again:
```
row.names(alpha_met_df) = alpha_met_df$Row.names
```

Get rid of the Row.names column:
```
alpha_met_df$Row.names = NULL
```
Check the output:
```
head(alpha_met_df)

            Shannon   Simpson InvSimpson    site_name    sex  race consent_age diagnosis
CSM5FZ3N_P 1.741766 0.6833787   3.158347 Cedars-Sinai Female White          43        CD
CSM5FZ3T_P 1.634631 0.7566009   4.108479 Cedars-Sinai Female White          76        CD
CSM5FZ4A_P 2.117305 0.8070880   5.183709 Cedars-Sinai Female White          47        UC
CSM5MCTZ_P 2.665711 0.8731621   7.884076 Cedars-Sinai   Male White          32        UC
CSM5MCU4_P 1.119533 0.3974255   1.659546 Cedars-Sinai Female White          53        CD
CSM5MCV1_P 1.478073 0.6609943   2.949803 Cedars-Sinai Female White          20        CD
```


#### Shannon alpha diversity: Univariate

Can do this with a for loop:
```
for (col in names(alpha_met_df)[-c(2:3)]){
  alpha_shannon_univ = anova(lm(as.formula(paste("Shannon ~ ", col)), data = alpha_met_df[-c(2:3)])) 
  print(alpha_shannon_univ)
}

Analysis of Variance Table

Response: Shannon
          Df Sum Sq Mean Sq F value Pr(>F)
Residuals 95 28.521 0.30022               
Analysis of Variance Table

Response: Shannon
          Df  Sum Sq Mean Sq F value Pr(>F)
site_name  4  1.6732 0.41831  1.4179 0.2344
Residuals 91 26.8477 0.29503               
Analysis of Variance Table

Response: Shannon
          Df  Sum Sq Mean Sq F value  Pr(>F)  
sex        1  0.9209 0.92088  3.1363 0.07981 .
Residuals 94 27.6001 0.29362                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
Analysis of Variance Table

Response: Shannon
          Df  Sum Sq Mean Sq F value Pr(>F)
race       3  0.9266 0.30886  1.0297 0.3832
Residuals 92 27.5944 0.29994               
Analysis of Variance Table

Response: Shannon
            Df  Sum Sq  Mean Sq F value Pr(>F)
consent_age  1  0.0105 0.010538  0.0347 0.8525
Residuals   94 28.5104 0.303302               
Analysis of Variance Table

Response: Shannon
          Df  Sum Sq Mean Sq F value Pr(>F)
diagnosis  2  0.4046 0.20231  0.6692 0.5146
Residuals 93 28.1163 0.30233               
```

* Question: Everyone likely knows what the Pr(>F) value means (the p-value), but what does the Df value tell us?


Alternatively, can do it one-by-one:
```
shannon_site = anova(lm(Shannon ~ site_name, data = alpha_met_df))
shannon_site

Analysis of Variance Table

Response: Shannon
          Df  Sum Sq Mean Sq F value Pr(>F)
site_name  4  1.6732 0.41831  1.4179 0.2344
Residuals 91 26.8477 0.29503  
```

* Question: Why did we need to use `alpha_met_df[-c(2:3)]` instead of the entire `alpha_met_df` in the models performed with a for loop, but not here in the writing out one variable at a time models?

```
shannon_sex = anova(lm(Shannon ~ sex, data = alpha_met_df))
shannon_sex

Analysis of Variance Table

Response: Shannon
          Df  Sum Sq Mean Sq F value  Pr(>F)  
sex        1  0.9209 0.92088  3.1363 0.07981 .
Residuals 94 27.6001 0.29362                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
```
shannon_race = anova(lm(Shannon ~ race, data = alpha_met_df))
shannon_race

Analysis of Variance Table

Response: Shannon
          Df  Sum Sq Mean Sq F value Pr(>F)
race       3  0.9266 0.30886  1.0297 0.3832
Residuals 92 27.5944 0.29994   
```
```
shannon_age = anova(lm(Shannon ~ consent_age, data = alpha_met_df))
shannon_age

Analysis of Variance Table

Response: Shannon
            Df  Sum Sq  Mean Sq F value Pr(>F)
consent_age  1  0.0105 0.010538  0.0347 0.8525
Residuals   94 28.5104 0.303302   
```
```
shannon_diagnosis = anova(lm(Shannon ~ diagnosis, data = alpha_met_df))
shannon_diagnosis

Analysis of Variance Table

Response: Shannon
          Df  Sum Sq Mean Sq F value Pr(>F)
diagnosis  2  0.4046 0.20231  0.6692 0.5146
Residuals 93 28.1163 0.30233    
```

* Question: What are univariate and multivariate tests and what benefit do both provide?

#### Shannon alpha diversity: Multivariate

Can do it without being verbose:
```
shannon_multi = anova(lm(Shannon ~ ., data = alpha_met_df[-c(2:3)]))
shannon_multi

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

Alternatively, can write out each variable in the model:
```
shannon_multi = anova(lm(Shannon ~ site_name + sex + race + consent_age + diagnosis, data = alpha_met_df))
shannon_multi

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

* Exercise: Now run univariate models for diagnosis on the Simpson and Inverse Simpson metrics.
    * Question: How might one easily do this with a shortcut in RStudio?


Now to look at all 3 measures together for diagnosis (univariate):
```
rbind(shannon_diagnosis, simpson_diagnosis, invsimpson_diagnosis)

Analysis of Variance Table

Response: Shannon
           Df  Sum Sq Mean Sq F value Pr(>F)
diagnosis   2    0.40  0.2023  0.6692 0.5146
Residuals  93   28.12  0.3023               
diagnosis1  2    0.01  0.0068  0.3803 0.6847
Residuals1 93    1.67  0.0180               
diagnosis2  2   12.64  6.3207  0.5316 0.5894
Residuals2 93 1105.77 11.8900   
```


### Beta Diversity: PERMANOVA tests using Bray-Curtis dissimilarity

Calculate Bray-Curtis dissimilarity:
```
bray_species = vegdist(species_filt, method = "bray")
```

* Question: What makes Bray-Curtis dissimilarity a good metric to use for this type of feature data?


#### Univariate

Can do this with a for loop:
```
for (col in names(metadata)){
  adonis_univ <- adonis(as.formula(paste("bray_species ~ ", col)), data = metadata)
  print(adonis_univ)
}

Call:
adonis(formula = as.formula(paste("bray_species ~ ", col)), data = metadata) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
site_name  4    1.6651 0.41628  1.3922 0.05767  0.016 *
Residuals 91   27.2095 0.29901         0.94233         
Total     95   28.8747                 1.00000         
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Call:
adonis(formula = as.formula(paste("bray_species ~ ", col)), data = metadata) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
sex        1    0.5925 0.59254  1.9694 0.02052  0.012 *
Residuals 94   28.2821 0.30087         0.97948         
Total     95   28.8747                 1.00000         
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Call:
adonis(formula = as.formula(paste("bray_species ~ ", col)), data = metadata) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
race       3    0.7778 0.25926 0.84891 0.02694  0.793
Residuals 92   28.0969 0.30540         0.97306       
Total     95   28.8747                 1.00000       

Call:
adonis(formula = as.formula(paste("bray_species ~ ", col)), data = metadata) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
consent_age  1    0.4683 0.46832  1.5497 0.01622  0.063 .
Residuals   94   28.4063 0.30219         0.98378         
Total       95   28.8747                 1.00000         
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Call:
adonis(formula = as.formula(paste("bray_species ~ ", col)), data = metadata) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
diagnosis  2    0.9399 0.46996  1.5646 0.03255  0.019 *
Residuals 93   27.9347 0.30037         0.96745         
Total     95   28.8747                 1.00000         
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

* Question: What does the R2 value tell us?

Alternatively, can do it one-by-one:
```
adonis_site_tax = adonis(bray_species ~ site_name, data = metadata)
adonis_site_tax

Call:
adonis(formula = bray_species ~ site_name, data = metadata) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
site_name  4    1.6651 0.41628  1.3922 0.05767  0.017 *
Residuals 91   27.2095 0.29901         0.94233         
Total     95   28.8747                 1.00000         
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
```
adonis_sex_tax = adonis(bray_species ~ sex, data = metadata)
adonis_sex_tax

Call:
adonis(formula = bray_species ~ sex, data = metadata) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
sex        1    0.5925 0.59254  1.9694 0.02052  0.009 **
Residuals 94   28.2821 0.30087         0.97948          
Total     95   28.8747                 1.00000          
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
```
adonis_age_tax = adonis(bray_species ~ consent_age, data = metadata)
adonis_age_tax

Call:
adonis(formula = bray_species ~ consent_age, data = metadata) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
consent_age  1    0.4683 0.46832  1.5497 0.01622   0.06 .
Residuals   94   28.4063 0.30219         0.98378         
Total       95   28.8747                 1.00000         
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
```
adonis_race_tax = adonis(bray_species ~ race, data = metadata)
adonis_race_tax

Call:
adonis(formula = bray_species ~ race, data = metadata) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
race       3    0.7778 0.25926 0.84891 0.02694  0.808
Residuals 92   28.0969 0.30540         0.97306       
Total     95   28.8747                 1.00000 
```
```
adonis_diagnosis_tax = adonis(bray_species ~ diagnosis, data = metadata)
adonis_diagnosis_tax

Call:
adonis(formula = bray_species ~ diagnosis, data = metadata) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
diagnosis  2    0.9399 0.46996  1.5646 0.03255  0.023 *
Residuals 93   27.9347 0.30037         0.96745         
Total     95   28.8747                 1.00000         
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

#### Multivariate

Can do it without being verbose:
```
adonis_multi_tax = adonis(bray_species ~ ., data = metadata)
adonis_multi_tax

Call:
adonis(formula = bray_species ~ ., data = metadata) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
site_name    4    1.6651 0.41628  1.4195 0.05767  0.014 *
sex          1    0.4848 0.48479  1.6531 0.01679  0.047 *
race         3    0.8897 0.29658  1.0113 0.03081  0.433  
consent_age  1    0.3545 0.35448  1.2088 0.01228  0.245  
diagnosis    2    0.8470 0.42349  1.4441 0.02933  0.030 *
Residuals   84   24.6335 0.29326         0.85312         
Total       95   28.8747                 1.00000         
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

Alternatively, can write out each variable in the model:
```
adonis_multi_tax = adonis(bray_species ~ site_name + sex + race + consent_age + diagnosis, data = metadata)
adonis_multi_tax

Call:
adonis(formula = bray_species ~ site_name + sex + race + consent_age + diagnosis, data = metadata) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
site_name    4    1.6651 0.41628  1.4195 0.05767  0.015 *
sex          1    0.4848 0.48479  1.6531 0.01679  0.051 .
race         3    0.8897 0.29658  1.0113 0.03081  0.452  
consent_age  1    0.3545 0.35448  1.2088 0.01228  0.227  
diagnosis    2    0.8470 0.42349  1.4441 0.02933  0.041 *
Residuals   84   24.6335 0.29326         0.85312         
Total       95   28.8747                 1.00000         
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```


## MGX pathways

### Feature table and metadata table creation and formatting

Download the mgx pathway realtive abundance data and put it into the data directory:
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
row.names(dna_path_unstratified)

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

   CSM67UH7    CSM79HHW    HSM67VDT    HSM6XRQB    HSM6XRR3    HSM7J4LP 
  1.0000000   1.0000005   1.0000002   0.9999997   1.0000002   0.9999997 
...
```

Filter for beta diversity:

Check the dimensions of dna pathways pre-filtering:
```
dim(dna_path_unstratified)

[1] 466  96
```
Filter:
```
dna_path_unstratified_filt = dna_path_unstratified[apply(dna_path_unstratified, 1, function(x) sum(x > 0.0001) > 0.1 * ncol(dna_path_unstratified)), ]
```
Check the dimensions of dna pathways post-filtering:
```
dim(dna_path_unstratified_filt)

[1] 309  96
```

Let's transpose the dataframes for easier use downstream, making the rows be the samples just like in the metadata:
```
dna_path_unstratified_filt = data.frame(t(dna_path_unstratified_filt), check.names = F)
dna_path_unstratified = data.frame(t(dna_path_unstratified), check.names = F)
```


### Beta Diversity: PERMANOVA tests using Bray-Curtis Dissimilarity

Calculate Bray-Curtis dissimilarity:
```
bray_dna_path_unstratified = vegdist(dna_path_unstratified_filt, method = "bray")
```

#### Univariate

Can do this with a for loop:
```
for (col in names(metadata)){
  adonis_univ <- adonis(as.formula(paste("bray_dna_path_unstratified ~ ", col)), data = metadata)
  print(adonis_univ)
}

Call:
adonis(formula = as.formula(paste("bray_dna_path_unstratified ~ ",      col)), data = metadata) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
site_name  4    0.1763 0.044076  1.1084 0.04646  0.302
Residuals 91    3.6187 0.039765         0.95354       
Total     95    3.7950                  1.00000       

Call:
adonis(formula = as.formula(paste("bray_dna_path_unstratified ~ ",      col)), data = metadata) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs  MeanSqs F.Model     R2 Pr(>F)
sex        1    0.0383 0.038313 0.95867 0.0101  0.441
Residuals 94    3.7566 0.039964         0.9899       
Total     95    3.7950                  1.0000       

Call:
adonis(formula = as.formula(paste("bray_dna_path_unstratified ~ ",      col)), data = metadata) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
race       3    0.1211 0.040381  1.0112 0.03192  0.402
Residuals 92    3.6738 0.039933         0.96808       
Total     95    3.7950                  1.00000       

Call:
adonis(formula = as.formula(paste("bray_dna_path_unstratified ~ ",      col)), data = metadata) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
consent_age  1    0.0367 0.036663   0.917 0.00966  0.449
Residuals   94    3.7583 0.039982         0.99034       
Total       95    3.7950                  1.00000       

Call:
adonis(formula = as.formula(paste("bray_dna_path_unstratified ~ ",      col)), data = metadata) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
diagnosis  2    0.1327 0.066335  1.6845 0.03496  0.062 .
Residuals 93    3.6623 0.039379         0.96504         
Total     95    3.7950                  1.00000         
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

Alternatively, can do it one-by-one:
```
adonis_site_dna_path = adonis(bray_dna_path_unstratified ~ site_name, data = metadata)
adonis_site_dna_path

Call:
adonis(formula = bray_dna_path_unstratified ~ site_name, data = metadata) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
site_name  4    0.1763 0.044076  1.1084 0.04646  0.294
Residuals 91    3.6187 0.039765         0.95354       
Total     95    3.7950                  1.00000       
```
```
adonis_sex_dna_path = adonis(bray_dna_path_unstratified ~ sex, data = metadata)
adonis_sex_dna_path

Call:
adonis(formula = bray_dna_path_unstratified ~ sex, data = metadata) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs  MeanSqs F.Model     R2 Pr(>F)
sex        1    0.0383 0.038313 0.95867 0.0101  0.437
Residuals 94    3.7566 0.039964         0.9899       
Total     95    3.7950                  1.0000       
```
```
adonis_age_dna_path = adonis(bray_dna_path_unstratified ~ consent_age, data = metadata)
adonis_age_dna_path

Call:
adonis(formula = bray_dna_path_unstratified ~ consent_age, data = metadata) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
consent_age  1    0.0367 0.036663   0.917 0.00966  0.453
Residuals   94    3.7583 0.039982         0.99034       
Total       95    3.7950                  1.00000       
```
```
adonis_race_dna_path = adonis(bray_dna_path_unstratified ~ race, data = metadata)
adonis_race_dna_path

Call:
adonis(formula = bray_dna_path_unstratified ~ race, data = metadata) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
race       3    0.1211 0.040381  1.0112 0.03192  0.421
Residuals 92    3.6738 0.039933         0.96808       
Total     95    3.7950                  1.00000       
```
```
adonis_diagnosis_dna_path = adonis(bray_dna_path_unstratified ~ diagnosis, data = metadata)
adonis_diagnosis_dna_path

Call:
adonis(formula = bray_dna_path_unstratified ~ diagnosis, data = metadata) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
diagnosis  2    0.1327 0.066335  1.6845 0.03496  0.054 .
Residuals 93    3.6623 0.039379         0.96504         
Total     95    3.7950                  1.00000         
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```


#### Multivariate

Can do it without being verbose:
```
adonis_multi_dna_path = adonis(bray_dna_path_unstratified ~ ., data = metadata)
adonis_multi_dna_path

Call:
adonis(formula = bray_dna_path_unstratified ~ ., data = metadata) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
site_name    4    0.1763 0.044076  1.1105 0.04646  0.296
sex          1    0.0311 0.031101  0.7836 0.00820  0.578
race         3    0.1243 0.041425  1.0437 0.03275  0.379
consent_age  1    0.0293 0.029283  0.7378 0.00772  0.636
diagnosis    2    0.1001 0.050034  1.2606 0.02637  0.205
Residuals   84    3.3339 0.039690         0.87852       
Total       95    3.7950                  1.00000       
```


* Question: Wait a minute! Diagnosis was near significant in the univariate model, but is now quite non-significant in the multivariate model. What happened?


Alternatively, can write out each variable in the model:
```
adonis_multi_dna_path = adonis(bray_dna_path_unstratified ~ site_name + sex + race + consent_age + diagnosis, data = metadata)
adonis_multi_dna_path

Call:
adonis(formula = bray_dna_path_unstratified ~ site_name + sex + race + consent_age + diagnosis, data = metadata) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
site_name    4    0.1763 0.044076  1.1105 0.04646  0.294
sex          1    0.0311 0.031101  0.7836 0.00820  0.580
race         3    0.1243 0.041425  1.0437 0.03275  0.409
consent_age  1    0.0293 0.029283  0.7378 0.00772  0.636
diagnosis    2    0.1001 0.050034  1.2606 0.02637  0.194
Residuals   84    3.3339 0.039690         0.87852       
Total       95    3.7950                  1.00000       
```

* Exercise: Run the same model a few times and compare the results.
    * Question: Are you seeing the same results across runs? If not, why is this?

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

[1] 96  5
```

Subset metadata and assign it to a new dataframe:
```
metadata_rna = subset(metadata, row.names(metadata) %in% names(rna_path))
```

Check the dimensions of the subsetted metadata:
```
dim(metadata_rna)

[1] 28  5
```

Check the output:
```
metadata_rna

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
row.names(rna_path_unstratified)

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

   HSM6XRQB    HSM7J4LP    MSM5LLDI    MSM6J2Q1    MSM79H8D    MSM9VZL5 
  1.0000001   0.9999992   0.9999997   1.0000000   1.0000001   0.9999998 
...
```

Filter for beta diversity:

Check the dimensions of rna pathways pre-filtering:
```
dim(rna_path_unstratified)

[1] 430  28
```
Filter:
```
rna_path_unstratified_filt = rna_path_unstratified[apply(rna_path_unstratified, 1, function(x) sum(x > 0.0001) > 0.1 * ncol(rna_path_unstratified)), ]
```
Check the dimensions of rna pathways post-filtering:
```
dim(rna_path_unstratified_filt)

[1] 234  28
```

Let's transpose the dataframes for easier use downstream, making the rows be the samples just like in the metadata:
```
rna_path_unstratified_filt = data.frame(t(rna_path_unstratified_filt), check.names = F)
rna_path_unstratified = data.frame(t(rna_path_unstratified), check.names = F)
```


### Beta Diversity: PERMANOVA tests using Bray-Curtis Dissimilarity

Calculate Bray-Curtis dissimilarity:
```
bray_rna_path_unstratified = vegdist(rna_path_unstratified_filt, method = "bray")
```

#### Univariate

Can do this with a for loop:
```
for (col in names(metadata)){
  adonis_univ <- adonis(as.formula(paste("bray_rna_path_unstratified ~ ", col)), data = metadata)
  print(adonis_univ)
}

Error in G * t(hat) : non-conformable arrays
```


* Question: Why does this result in an error? 
    * Hint: Check the formula call and remember what we subsetted for this rna data.


Alternatively, can do it one-by-one:
```
adonis_site_rna_path = adonis(bray_rna_path_unstratified ~ site_name, data = metadata_rna)
adonis_site_rna_path

Call:
adonis(formula = bray_rna_path_unstratified ~ site_name, data = metadata_rna) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
site_name  4   0.33933 0.084834  1.0482 0.15419  0.385
Residuals 23   1.86136 0.080929         0.84581       
Total     27   2.20069                  1.00000   
```
```
adonis_sex_rna_path = adonis(bray_rna_path_unstratified ~ sex, data = metadata_rna)
adonis_sex_rna_path

Call:
adonis(formula = bray_rna_path_unstratified ~ sex, data = metadata_rna) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
sex        1   0.12339 0.123389  1.5444 0.05607  0.062 .
Residuals 26   2.07730 0.079896         0.94393         
Total     27   2.20069                  1.00000         
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1   
```
```
adonis_age_rna_path = adonis(bray_rna_path_unstratified ~ consent_age, data = metadata_rna)
adonis_age_rna_path

Call:
adonis(formula = bray_rna_path_unstratified ~ consent_age, data = metadata_rna) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
consent_age  1   0.08308 0.083081  1.0201 0.03775  0.295
Residuals   26   2.11761 0.081447         0.96225       
Total       27   2.20069                  1.00000       
```
```
adonis_race_rna_path = adonis(bray_rna_path_unstratified ~ race, data = metadata_rna)
adonis_race_rna_path

 Call:
adonis(formula = bray_rna_path_unstratified ~ race, data = metadata_rna) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
race       3   0.13184 0.043945 0.50979 0.05991  0.901
Residuals 24   2.06886 0.086202         0.94009       
Total     27   2.20069                  1.00000       
```
```
adonis_diagnosis_rna_path = adonis(bray_rna_path_unstratified ~ diagnosis, data = metadata_rna)
adonis_diagnosis_rna_path

Call:
adonis(formula = bray_rna_path_unstratified ~ diagnosis, data = metadata_rna) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
diagnosis  2   0.09497 0.047483 0.56373 0.04315  0.988
Residuals 25   2.10573 0.084229         0.95685       
Total     27   2.20069                  1.00000    
```

Solution to the for loop question earlier:
```
for (col in names(metadata_rna)){
  adonis_univ <- adonis(as.formula(paste("bray_rna_path_unstratified ~ ", col)), data = metadata_rna)
  print(adonis_univ)
}

Call:
adonis(formula = as.formula(paste("bray_rna_path_unstratified ~ ",      col)), data = metadata_rna) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
site_name  4   0.33933 0.084834  1.0482 0.15419  0.366
Residuals 23   1.86136 0.080929         0.84581       
Total     27   2.20069                  1.00000       

Call:
adonis(formula = as.formula(paste("bray_rna_path_unstratified ~ ",      col)), data = metadata_rna) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
sex        1   0.12339 0.123389  1.5444 0.05607  0.049 *
Residuals 26   2.07730 0.079896         0.94393         
Total     27   2.20069                  1.00000         
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Call:
adonis(formula = as.formula(paste("bray_rna_path_unstratified ~ ",      col)), data = metadata_rna) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
race       3   0.13184 0.043945 0.50979 0.05991  0.894
Residuals 24   2.06886 0.086202         0.94009       
Total     27   2.20069                  1.00000       

Call:
adonis(formula = as.formula(paste("bray_rna_path_unstratified ~ ",      col)), data = metadata_rna) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
consent_age  1   0.08308 0.083081  1.0201 0.03775  0.331
Residuals   26   2.11761 0.081447         0.96225       
Total       27   2.20069                  1.00000       

Call:
adonis(formula = as.formula(paste("bray_rna_path_unstratified ~ ",      col)), data = metadata_rna) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
diagnosis  2   0.09497 0.047483 0.56373 0.04315  0.989
Residuals 25   2.10573 0.084229         0.95685       
Total     27   2.20069                  1.00000       
```

#### Multivariate

Can do it without being verbose:
```
adonis_multi_rna_path = adonis(bray_rna_path_unstratified ~ ., data = metadata_rna)
adonis_multi_rna_path

Call:
adonis(formula = bray_rna_path_unstratified ~ ., data = metadata_rna) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
site_name    4   0.33933 0.084834 0.95180 0.15419  0.487
sex          1   0.09628 0.096283 1.08026 0.04375  0.366
race         3   0.15885 0.052948 0.59406 0.07218  0.806
consent_age  1   0.06621 0.066206 0.74281 0.03008  0.588
diagnosis    2   0.11395 0.056976 0.63925 0.05178  0.940
Residuals   16   1.42607 0.089129         0.64801       
Total       27   2.20069                  1.00000       
```

Alternatively, can write out each variable in the model:
```
adonis_multi_rna_path = adonis(bray_rna_path_unstratified ~ site_name + sex + race + consent_age + diagnosis, data = metadata_rna)
adonis_multi_rna_path

Call:
adonis(formula = bray_rna_path_unstratified ~ site_name + sex +      race + consent_age + diagnosis, data = metadata_rna) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
site_name    4   0.33933 0.084834 0.95180 0.15419  0.499
sex          1   0.09628 0.096283 1.08026 0.04375  0.393
race         3   0.15885 0.052948 0.59406 0.07218  0.836
consent_age  1   0.06621 0.066206 0.74281 0.03008  0.585
diagnosis    2   0.11395 0.056976 0.63925 0.05178  0.933
Residuals   16   1.42607 0.089129         0.64801       
Total       27   2.20069                  1.00000       
```



## RNA/DNA pathway ratios

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
```
dim(metadata)
```

Check the dimensions of `metadata` pre-subsetting:
```
dim(metadata)

[1] 96  5
```

Subset metadata and assign it to a new dataframe:
```
metadata_rna_dna = subset(metadata, row.names(metadata) %in% names(rna_dna_path))
```

Check the dimensions of the subsetted metadata:
```
dim(metadata_rna_dna)

[1] 12  5
```

Check the output:
```
metadata_rna_dna

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
row.names(rna_dna_path_unstratified)

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

[1] 472  12
```
Filter:
```
rna_dna_path_unstratified_filt = subset(rna_dna_path_unstratified, row.names(rna_dna_path_unstratified) %in% names(dna_path_unstratified_filt))
```
Check the dimensions of RNA/DNA pathways post-filtering:
```
dim(rna_dna_path_unstratified_filt)

[1] 309  12
```
Check the dimensions of dna pathways post-filtering to make sure they line up:
```
dim(dna_path_unstratified_filt)

[1]  96 309
```
Yes, they have the same number of pathways.


Let's transpose the dataframes for easier use downstream, making the rows be the samples just like in the metadata:
```
rna_dna_path_unstratified_filt = data.frame(t(rna_dna_path_unstratified_filt), check.names = F)
rna_dna_path_unstratified = data.frame(t(rna_dna_path_unstratified), check.names = F)
```

log transform the RNA/DNA ratio:
```
rna_dna_path_unstratified_filt_log = log2(rna_dna_path_unstratified_filt + 1)
```


### Beta Diversity: PERMANOVA tests using Euclidean distances

Calculate Euclidean distances:
```
euclidean_rna_dna_path_unstratified = vegdist(rna_dna_path_unstratified_filt_log, method = "euclidean")
```

* Question: Why are we using Euclidean distances on this dataset as opposed to Bray-Curtis dissimilarities that we used on all of the other feature tables?


#### Univariate

Can do this with a for loop:
```
for (col in names(metadata_rna_dna)){
  adonis_univ <- adonis(as.formula(paste("euclidean_rna_dna_path_unstratified ~ ", col)), data = metadata_rna_dna)
  print(adonis_univ)
}

Call:
adonis(formula = as.formula(paste("euclidean_rna_dna_path_unstratified ~ ",      col)), data = metadata_rna_dna) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)
site_name  2    257.63  128.82  1.0093 0.1832  0.418
Residuals  9   1148.66  127.63         0.8168       
Total     11   1406.29                 1.0000       

Call:
adonis(formula = as.formula(paste("euclidean_rna_dna_path_unstratified ~ ",      col)), data = metadata_rna_dna) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
sex        1    148.46  148.46  1.1803 0.10557  0.204
Residuals 10   1257.84  125.78         0.89443       
Total     11   1406.29                 1.00000       

Error in `contrasts<-`(`*tmp*`, value = contr.funs[1 + isOF[nn]]) : 
  contrasts can be applied only to factors with 2 or more levels
```


* Question: We got an error, what could be the cause?
    * Hint: Let's investigate the metadata for this set of samples again.


```
metadata_rna_dna

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

We only have one race now in this dataset, so let's remove that column from the metadata in the formula call. This wouldn't be a problem in the one-by-one method below, but is in the loop that goes through each column one at a time. Rememeber how we did this before for the taxonomy data in the alpha diversity models?
```
for (col in names(metadata_rna_dna)[-3]){
  adonis_univ <- adonis(as.formula(paste("euclidean_rna_dna_path_unstratified ~ ", col)), data = metadata_rna_dna[-3])
  print(adonis_univ)
}

Call:
adonis(formula = as.formula(paste("euclidean_rna_dna_path_unstratified ~ ",      col)), data = metadata_rna_dna[-3]) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)
site_name  2    257.63  128.82  1.0093 0.1832  0.424
Residuals  9   1148.66  127.63         0.8168       
Total     11   1406.29                 1.0000       

Call:
adonis(formula = as.formula(paste("euclidean_rna_dna_path_unstratified ~ ",      col)), data = metadata_rna_dna[-3]) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
sex        1    148.46  148.46  1.1803 0.10557  0.208
Residuals 10   1257.84  125.78         0.89443       
Total     11   1406.29                 1.00000       

Call:
adonis(formula = as.formula(paste("euclidean_rna_dna_path_unstratified ~ ",      col)), data = metadata_rna_dna[-3]) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
consent_age  1    123.86  123.86 0.96584 0.08808  0.457
Residuals   10   1282.43  128.24         0.91192       
Total       11   1406.29                 1.00000       

Call:
adonis(formula = as.formula(paste("euclidean_rna_dna_path_unstratified ~ ",      col)), data = metadata_rna_dna[-3]) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
diagnosis  2    332.89  166.44  1.3956 0.23671  0.056 .
Residuals  9   1073.41  119.27         0.76329         
Total     11   1406.29                 1.00000         
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

Alternatively, can do it one-by-one:
```
adonis_site_rna_dna_path = adonis(euclidean_rna_dna_path_unstratified ~ site_name, data = metadata_rna_dna)
adonis_site_rna_dna_path

Call:
adonis(formula = euclidean_rna_dna_path_unstratified ~ site_name,      data = metadata_rna_dna) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)
site_name  2    257.63  128.82  1.0093 0.1832  0.394
Residuals  9   1148.66  127.63         0.8168       
Total     11   1406.29                 1.0000   
```
```
adonis_sex_rna_dna_path = adonis(euclidean_rna_dna_path_unstratified ~ sex, data = metadata_rna_dna)
adonis_sex_rna_dna_path

Call:
adonis(formula = euclidean_rna_dna_path_unstratified ~ sex, data = metadata_rna_dna) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
sex        1    148.46  148.46  1.1803 0.10557  0.206
Residuals 10   1257.84  125.78         0.89443       
Total     11   1406.29                 1.00000    
```
```
adonis_age_rna_dna_path = adonis(euclidean_rna_dna_path_unstratified ~ consent_age, data = metadata_rna_dna)
adonis_age_rna_dna_path

Call:
adonis(formula = euclidean_rna_dna_path_unstratified ~ consent_age,      data = metadata_rna_dna) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
consent_age  1    123.86  123.86 0.96584 0.08808  0.434
Residuals   10   1282.43  128.24         0.91192       
Total       11   1406.29                 1.00000    
```
```
adonis_diagnosis_rna_dna_path = adonis(euclidean_rna_dna_path_unstratified ~ diagnosis, data = metadata_rna_dna)
adonis_diagnosis_rna_dna_path

Call:
adonis(formula = euclidean_rna_dna_path_unstratified ~ diagnosis,      data = metadata_rna_dna) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
diagnosis  2    332.89  166.44  1.3956 0.23671   0.06 .
Residuals  9   1073.41  119.27         0.76329         
Total     11   1406.29                 1.00000         
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

#### Multivariate

Can do it without being verbose:
```
adonis_multi_rna_dna_path = adonis(euclidean_rna_dna_path_unstratified ~ ., data = metadata_rna_dna[-3])
adonis_multi_rna_dna_path

Call:
adonis(formula = euclidean_rna_dna_path_unstratified ~ ., data = metadata_rna_dna[-3]) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
site_name    2    257.63  128.82 1.00288 0.18320  0.446
sex          1    128.11  128.11 0.99736 0.09110  0.420
consent_age  1    112.95  112.95 0.87936 0.08032  0.616
diagnosis    2    265.36  132.68 1.03293 0.18869  0.407
Residuals    5    642.24  128.45         0.45669       
Total       11   1406.29                 1.00000       
```

Alternatively, can write out each variable in the model:
```
adonis_multi_rna_dna_path = adonis(euclidean_rna_dna_path_unstratified ~ site_name + sex + consent_age + diagnosis, data = metadata_rna_dna)
adonis_multi_rna_dna_path

Call:
adonis(formula = euclidean_rna_dna_path_unstratified ~ site_name +      sex + consent_age + diagnosis, data = metadata_rna_dna) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
site_name    2    257.63  128.82 1.00288 0.18320  0.462
sex          1    128.11  128.11 0.99736 0.09110  0.436
consent_age  1    112.95  112.95 0.87936 0.08032  0.604
diagnosis    2    265.36  132.68 1.03293 0.18869  0.432
Residuals    5    642.24  128.45         0.45669       
Total       11   1406.29                 1.00000       
```

Diagnosis goes from being significant in the univariate model to non-significant in the multivariate.


Now let's look at all of the PERMANOVA results together:
```
adonis_multi_tax

Call:
adonis(formula = bray_species ~ site_name + sex + race + consent_age +      diagnosis, data = metadata) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
site_name    4    1.6651 0.41628  1.4195 0.05767  0.015 *
sex          1    0.4848 0.48479  1.6531 0.01679  0.051 .
race         3    0.8897 0.29658  1.0113 0.03081  0.452  
consent_age  1    0.3545 0.35448  1.2088 0.01228  0.227  
diagnosis    2    0.8470 0.42349  1.4441 0.02933  0.041 *
Residuals   84   24.6335 0.29326         0.85312         
Total       95   28.8747                 1.00000         
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
```
adonis_multi_dna_path

Call:
adonis(formula = bray_dna_path_unstratified ~ site_name + sex +      race + consent_age + diagnosis, data = metadata) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
site_name    4    0.1763 0.044076  1.1105 0.04646  0.294
sex          1    0.0311 0.031101  0.7836 0.00820  0.580
race         3    0.1243 0.041425  1.0437 0.03275  0.409
consent_age  1    0.0293 0.029283  0.7378 0.00772  0.636
diagnosis    2    0.1001 0.050034  1.2606 0.02637  0.194
Residuals   84    3.3339 0.039690         0.87852       
Total       95    3.7950                  1.00000       
```
```
adonis_multi_rna_path

Call:
adonis(formula = bray_rna_path_unstratified ~ site_name + sex +      race + consent_age + diagnosis, data = metadata_rna) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
site_name    4   0.33933 0.084834 0.95180 0.15419  0.499
sex          1   0.09628 0.096283 1.08026 0.04375  0.393
race         3   0.15885 0.052948 0.59406 0.07218  0.836
consent_age  1   0.06621 0.066206 0.74281 0.03008  0.585
diagnosis    2   0.11395 0.056976 0.63925 0.05178  0.933
Residuals   16   1.42607 0.089129         0.64801       
Total       27   2.20069                  1.00000       
```
```
adonis_multi_rna_dna_path

Call:
adonis(formula = euclidean_rna_dna_path_unstratified ~ site_name +      sex + consent_age + diagnosis, data = metadata_rna_dna) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
site_name    2    257.63  128.82 1.00288 0.18320  0.462
sex          1    128.11  128.11 0.99736 0.09110  0.436
consent_age  1    112.95  112.95 0.87936 0.08032  0.604
diagnosis    2    265.36  132.68 1.03293 0.18869  0.432
Residuals    5    642.24  128.45         0.45669       
Total       11   1406.29                 1.00000       
```

* Question: What are some possible reasons the taxonomy data gives the only significant results?


# Mantel tests

To be continued...