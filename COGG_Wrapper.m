%*************************COGG and COGG-CCA********************************
%Author: Aritra Bose, PhD student,Computer Science, Purdue University
%Last Edit on Oct 13, 2016
%The following code implements COGG or Correlation Optimization of Genetics
%and Geodemographics. This is a wrapper code which processes the input
%files and prepares it into a vector U containing the principal component
%and a matrix G containing the Geodemographic Matrix. It then computes COGG
%and returns two vectors for PC1 and PC2 respectively. The vectors look like:
%[OriginalSquaredCorrelation RandomlyPermutedCorrelations COGGResult]
%this vector can be used to plot a histogram and see the statistical 
%significance of COGG. The last value in this vector is the result that COGG gives. 
%We also use this as a wrapper for COGG-CCA which is extending COGG to a 
%Canonical Correlation Analysis setup. CCA finds the linear combinations of 
%the columns of the Principal Component Matrix and Geodemographic Matrix 
%and return the maximum correlation with each other. COGG-CCA would return
%the same type of vector, with the first element being the Original Correlation
%without COGG-CCA, followed by Random Permutations of Caste and
%Languages and finally COGG-CCA's result showing statistical significance. 

%List of Inputs: 
% --> File containing Principal Components (we use EIGENSTRAT output)
% --> File containing external information along with geographical
%     coordinates. In our case we use Caste and Language information as
%     external information
% --> value of p, the top principal components to be considered. 
clear; 
close all; 
%Clearing all existing variables
PCs_PathToFile = 'NormIN_AllRep/NormIN_AllRep_targetpop_pca_toplot.txt'; 
G_PathToFile = 'Norm_IN_AllRep.xlsx';
%Defining the paths to the files containing the principal components and
%Geodemographics. The Geodemographic file that we use contain the sample
%IDs, corresponding Latitude and Longitude information in decimals, the
%Caste and Language group they belong to. 
PCs = table2cell(readtable(PCs_PathToFile,'Delimiter',' '));
[r,c] = size(PCs); 
G = table2cell(readtable(G_PathToFile));
PC1and2 = [cell2mat(PCs(:,2)) cell2mat(PCs(:,3))];
A = PCs(:,2:12);
p = 8; pPCs = [];
for i = 1:p
    pPCs = [pPCs cell2mat(A(:,i))];
end
GeoMat = [G(:,1) G(:,7) G(:,8) G(:,6) G(:,5)];
coords = zeros(length(PCs),2);
castes = cell(length(PCs),1);
languages = cell(length(PCs),1);
%Defining the empty vectors for caste, languages and geographical
%coordinates
for i = 1: length(GeoMat(:,1))
    idx = strmatch(GeoMat(i,1), PCs(:,13),'exact');
    coords(idx,1) = cell2mat(GeoMat(i,2)) ;
    coords(idx,2) = cell2mat(GeoMat(i,3));
    castes(idx) = GeoMat(i,4);
    languages(idx) = GeoMat(i,5);
end
%We matched the populations with the corresponding IDs in the two input
%files and for each sample found their language, caste and Geographical
%Coordinates.
UniqCastes = unique(castes);
NumCastes = length(UniqCastes);
UniqLanguages = unique(languages);
NumLangs = length(UniqLanguages);
langmat = []; castemat = [];
for i = 1:length(UniqCastes)
    tmpcaste = strcmp(UniqCastes(i),castes);
    castemat = [castemat tmpcaste];
end
Caste_Matrix = castemat;
for i = 1:length(UniqLanguages)
    tmplang = strcmp(UniqLanguages(i),languages);
    langmat = [langmat tmplang];
end
Language_Matrix = langmat; 
%in the above we encode the caste and language matrix as switch values,
%like 0 and 1 binary value for each sample belonging to a particular caste
%and language, we use the value 1 for that caste and language and 0
%everywhere else. The Language_Matrix and Caste_Matrix, matrices are sparse
%matrices. 
if(numel(find(coords == 0))> 0)
    fprintf('\nPopulation Missing! \n')
    unique(PCs(find(coords(:,1) == 0),13))
    return
end
%finding for missing values, for double-checking. 
%GeoDem = [coords Caste_Matrix Language_Matrix];
%We read the files and arrange the Principal Components in the matrix
%PC1and2 and the Geodemographic information in the GeoDem matrix. For,
%COGG-CCA we also define pPCs which contain the top "p" Principal
%Components. "p" is user defined. 
[COGG_PC1_Out, COGG_PC2_Out] = ComputeCOGG(PC1and2,coords,Caste_Matrix,Language_Matrix,NumCastes,NumLangs); 
figure(1);
histogram(COGG_PC1_Out); 
title('Distribution of correlation values showing statistical significance for PC1');xlabel('COGG values');ylabel('frequency'); 
figure(2);
histogram(COGG_PC2_Out); 
title('Distribution of correlation values showing statistical significance for PC2');xlabel('COGG values');ylabel('frequency');
%ComputeCOGG functions computes the maximum squared correlation between 
%PC1 OR PC2 and GeoDem. 
COGG_CCA_Out = ComputeCOGG_CCA(pPCs,coords,Caste_Matrix,Language_Matrix,NumCastes,NumLangs);
%ComputeCOGG functions computes the maximum squared correlation between 
%PC1 OR PC2 and GeoDem. 
figure(3);
histogram(COGG_CCA_Out); 
title('Distribution of correlation values showing statistical significance');xlabel('COGG-CCA values');ylabel('frequency');