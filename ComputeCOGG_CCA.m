function COGGCCA_Out = ComputeCOGG_CCA(pPCs,coords,CasteMat,LangMat,NumCastes,NumLangs)
%This function takes input the p Principal Components, Geographical
%Coordinates, Caste and Language information and computes the maximum
%Canonical Correlation among the Principal Component matrix and
%Geodemographic Matrix. In order to do that we solve for the coefficient
%values for which we get maximum correlation and plug it in and return the
%max correlation. We also check for statistical significance by randomly
%permuting rows of Castes and Languages in the Geodemographic Matrix. We do
%this 1000 times. 
ZPCs = zscore(pPCs);
ZG1 = zscore(coords(:,1));ZG2 = zscore(coords(:,2));
OrigCCACOGG = COGGCCAfun(ZPCs,[ZG1 ZG2]);
fprintf('Before COGG-CCA: %4.6f\n', OrigCCACOGG); 
%Computing the CCA between Principal Components and Geographical
%coordinates, before using Caste and Languages
for i = 1:NumCastes%Zscores for Caste Matrix
    Zcastes{i} = zscore(CasteMat(:,i))';  
end
for i = 1:NumLangs%Zscores for Language Matrix 
    Zlangs{i} = zscore(LangMat(:,i))';  
end
ZCompFac = [ZG1 ZG2 Zcastes{1,1}' Zcastes{1,2}' Zcastes{1,3}' Zlangs{1,1}' Zlangs{1,2}' Zlangs{1,3}' Zlangs{1,4}'];
CSC = COGGCCAfun(ZPCs, ZCompFac); %Computing COGG-CCA value with Principal Components and Geodemographic Matrix
CCACOGGcorr = CSC; 
fprintf('Fitting the solution for COGG-CCA: %4.6f\n', CCACOGGcorr);
RandSC = []; 
randcutoff = 1; 
for iter = 1:1000 %We do this a 1000 times
    [r,c] = size(LangMat);
    for i = 1:c%Randomly Permuting the rows of Language matrix
        randcutidx = randperm(floor(r*randcutoff));
        randidx = [randcutidx ((floor(r*randcutoff)+1):r) ];
        randlangmat(:,i) = LangMat(randidx,i);
    end
    for i = 1:NumLangs
       Zrandlangs{i} = zscore(randlangmat(:,i))';  
    end
    [r,c] = size(CasteMat);
    randcastevec = zeros(r,c);
    for i = 1:c%Randomly Permuting the rows of Caste matrix
        randcutidx = randperm(floor(r*randcutoff));
        randidx = [randcutidx ((floor(r*randcutoff)+1):r) ];
        randcastevec(:,i) = CasteMat(randidx,i);
    end
    for i = 1:NumCastes
       Zcastevecrand{i} = zscore(randcastevec(:,i))';  
    end
    ZRandCompFac = [ZG1 ZG2 Zcastevecrand{1,1}' Zcastevecrand{1,2}' Zcastevecrand{1,3}' Zrandlangs{1,1}' Zrandlangs{1,2}' Zrandlangs{1,3}' Zrandlangs{1,4}'];
    %ZRandCompFac is the matrix with geographical coordinates and randomly
    %permuted caste and language information
    CSC = COGGCCAfun(ZPCs,ZRandCompFac);
    %Saving each value of the resultant squared correlation in a vector after each iteration
    RandSC = [RandSC CSC];
end
COGGCCA_Out = [OrigCCACOGG RandSC CCACOGGcorr]; 