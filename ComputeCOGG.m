function [COGG_PC1_Out,COGG_PC2_Out] = ComputeCOGG(PC1and2,coords,CasteMat,LangMat,NumCastes,NumLangs)

origval = corr(coords,PC1and2);
%We calculate for original correlation with PC1and2 and geographical
%coordinates without encoding Caste and Language information.
sq_origval = origval.^2;
%We take the squared correlation as we want to capture the magnitude, not
%the direction of the correlation.
[maxval,maxidx] = max(sq_origval);
if maxval(1) > maxval(2)
    if maxidx == 1
        PC1_Lat = sq_origval(1,1);
        PC2_Long = sq_origval(2,2);
        OSC =[PC1_Lat PC2_Long];
    else
        PC1_Long = sq_origval(1,2);
        PC2_Lat = sq_origval(2,1);
        OSC = [PC1_Long PC2_Lat]; 
    end
else
    if maxidx == 1
        PC2_Lat = sq_origval(2,1); 
        PC1_Long = sq_origval(1,2);
        OSC =[PC1_Long PC2_Lat];
    else
        PC1_Lat = sq_origval(1,1);
        PC2_Long = sq_origval(2,2); 
        OSC = [PC1_Lat PC2_Long]; 
    end
end
fprintf('Before COGG, PC2: %4.6f\n', OSC(2) );
fprintf('Before COGG, PC1: %4.6f\n', OSC(1) );
%Here we check whether PC1 is related to Longitude and Latitude and same
%goes for PC2. We save the values to Original Squared Correlation variable,
%OSC.
%Now, we convert all of the encoded information in the GeoDem matrix to
%Z-scores for normalizing them, before using. 
ZPC1 = zscore(PC1and2(:,1));ZPC2 = zscore(PC1and2(:,2));
ZG1 = zscore(coords(:,1));ZG2 = zscore(coords(:,2));
for i = 1:NumCastes
    Zcastes{i} = zscore(CasteMat(:,i))';  
end
for i = 1:NumLangs
    Zlangs{i} = zscore(LangMat(:,i))';  
end
ZPCs = [ZPC1 ZPC2];
ZCompFac = [ZG1 ZG2 Zcastes{1,1}' Zcastes{1,2}' Zcastes{1,3}' Zlangs{1,1}' Zlangs{1,2}' Zlangs{1,3}' Zlangs{1,4}'];
%Here, we create the Geodemographic matrix with the Z-values of encoded
%information from Geographical Coordinates,Castes,Languages
%Now, we compute the original result for COGG
CSC = COGGfun(ZPCs,ZCompFac); 
COGGcorr = CSC;
fprintf('Fitting the solution for PC1: %4.6f\n', COGGcorr(1));
fprintf('Fitting the solution for PC2: %4.6f\n', COGGcorr(2));

%To check statistical significance now, we permute all Castes and Language
%encoded information and compute the correlation 
RandSC1 = [];RandSC2 = [];
randcutoff =1; %this value decides how much of the Castes and Languages to
%permute
for iter = 1:1000 %We do this a 1000 times
    [r,c] = size(LangMat);
    for i = 1:c
        randcutidx = randperm(floor(r*randcutoff));
        randidx = [randcutidx ((floor(r*randcutoff)+1):r) ];
        randlangmat(:,i) = LangMat(randidx,i);
    end
    for i = 1:NumLangs
       Zrandlangs{i} = zscore(randlangmat(:,i))';  
    end
    [r,c] = size(CasteMat);
    randcastevec = zeros(r,c);
    for i = 1:c
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
    CSC = COGGfun(ZPCs,ZRandCompFac);
    %Saving each value of the resultant squared correlation in a vector after each iteration
    RandSC1 = [RandSC1 CSC(1)];
    RandSC2 = [RandSC2 CSC(2)];
end
COGG_PC1_Out = [OSC(1) RandSC1 COGGcorr(2)];
COGG_PC2_Out = [OSC(2) RandSC2 COGGcorr(1)];