function CCACOGGval = COGGCCAfun(ZPCs,ZCompFac)
[~,c1] = size(ZPCs); 
[~,c2] = size(ZCompFac);
covmat = [];
%computes the covariance matrix, that is the covariance between each
%Principal Component vector and each column in the Geodemographic Matrix
for i = 1:c1 
    tmpcovmat = []; 
    for j = 1:c2
        tmpcov = cov(ZPCs(:,i),ZCompFac(:,j));
        tmpcovmat = [tmpcovmat tmpcov(1,2)];
    end
    covmat = [covmat;tmpcovmat];
end
%Computes the SVD of the principal component matrix
[U_xx2,S_xx2,V_xx2] = svd(cov(ZPCs));
eig_XX2 = diag(S_xx2); %getting the singular values out
neweig_XX2 = 1./eig_XX2;%reciprocal of each singular value
newSxx2 = diag(sqrt(neweig_XX2));
SigmaXX2 = U_xx2*newSxx2*V_xx2';%forming the new Princpal Component matrix's inverse

%Computes the SVD of the Geodemographic matrix
[U_yy,S_yy,V_yy] = svd(cov(ZCompFac)); 
eig_YY = diag(S_yy);%getting the singular values out
neweig_YY = 1./eig_YY; %reciprocal of each singular value
newSyy = diag(sqrt(neweig_YY));
SigmaYY = U_yy*newSyy*V_yy';%forming the inverse of the Geodemographic Matrix

BigSigma2 = SigmaYY*covmat'*SigmaXX2; %Forming the Sigma matrix, which is the solution to CCA

[ls,ss,rs] = svd(BigSigma2); %Computing the SVD of the Big Sigma matrix
d = diag(ss);
alpha = rs*SigmaXX2; %finding the solution for alphas
beta = SigmaYY*ls; %finding the solution for betas
ccacogg1 = 0; 
ccacogg2 = 0;
coeffPCs = [];
coeffeverything = [];
for i = 1:length(rs)
    sum = alpha(i,1)*ZPCs(:,i); %computing the COGG-CCA values for Principal Component Matrix
    coeffPCs = [coeffPCs sum]; 
    ccacogg1 = ccacogg1 + sum;
end
for i = 1:length(ls)
    sum = beta(i,1)*ZCompFac(:,i);%computing the COGG-CCA values for Geodemographic Matrix
    coeffeverything = [coeffeverything sum];
    ccacogg2 = ccacogg2 + sum;
end
CCACOGGval = corr(ccacogg1,ccacogg2)^2;%Computing squared correlation of COGG-CCA
