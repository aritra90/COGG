# COGG
#COGG
COGG or Correlation Optimization of Genetics and Geodemographics is a novel method where we find the 
maximum correlation of the principal components obtained from genotype data and the Geodemographic matrix,
which consists of Geographical coordinates and external information that influence the genetics of a particular 
region. In other words, if any external factors such as society, languages, occupations, natural phenomenon
has created genetic stratification in the population. 
In our case, we have chosen data from the Indian subcontinent. India is a land of diversity where the Genetics 
of the country has been largely influenced by prevalence of caste system and presence of different language 
families. 
COGG maximizes the squared correlation of the vector of Principal Component and Geodemographic matrix.


List of Inputs: 
--> File containing Principal Components (we use EIGENSTRAT output)
--> File containing external information along with geographical
    coordinates. In our case we use Caste and Language information as
    external information
--> value of p, the top principal components to be considered. 

#COGG-CCA
This also computes a Canonical Correlation Analysis of CCA, which we name COGG-CCA. COGG is one-sided CCA, 
we extend both sides in CCA. We now consider top "p" principal components, instead of a vector for the 
Principal Components. We want to maximize the squared correlation of the Principal Component matrix and 
the Geodemographic matrix. Canonical Correlation Analysis finds linear combinations of the variables 
in the Genetics and Geodemographic matrices.

Running COGG: 

You have to clone/download the codes and run the COGG_Wrapper.m 
Edit the COGG_Wrapper.m to include the paths to your files containing the Principal Components in the following
format: SampleID PC1 PC2 PC3 PC4 PC5 .. PCp Population
Where, SampleID is the id of each sample under study and the last column contains the population to which the 
samples belong. The PC1-PCp, are the number of principal components under consideration. 

Also include the link to the file containing external information of geographical coordinates and external 
information which influence the genetics of those populations. Our methods needs a one-to-one relationship
with the external factors and the sample. That is, each sample can belong to only one language group and caste 
group in our case. 

Also you have to edit the value of p and enter the desired number of principal components. 

Run it as COGG_Wrapper from the command line. 
