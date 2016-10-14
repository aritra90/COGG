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