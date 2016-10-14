function CSC = COGGfun(ZPCs,ZCompFac)
cov_compFac = cov(ZCompFac);
%we calculate the covariance of Geodemographic matrix
[~,c1] = size(ZPCs);
[~,c2] = size(ZCompFac);
covmat = []; 
for i = 1:c1 
    tmpcovmat = []; 
    for j = 1:c2
        tmpcov = cov(ZPCs(:,i),ZCompFac(:,j));
        tmpcovmat = [tmpcovmat tmpcov(1,2)];
    end
    covmat = [covmat;tmpcovmat];
end
%Computing the SVD of the covariance matrix
[U1,S1,V1] = svd(cov_compFac); 
el1 = diag(S1); %Vector of singular values
idx = find(el1 < 0.05);%Check whether the values are close to 0 leading 
%to poorly conditioned matrix during the inverse. We only keep values which
%are not extremely close to 0.
if idx > 0
    cut = idx-1;
else
    cut = c2; %If all the values are not close to 0, we consider all the components
end
newdiag1 = 1./el1(1:cut);%reciprocal of the singular values
newS1 = diag(newdiag1);
newD1 = U1(:,1:cut)*newS1*V1(:,1:cut)'; %The inverse matrix
cogg1 = 0; cogg2 = 0;
sol_coeffs1 = covmat(1,:)'\newD1;
sol_coeffs2 = covmat(2,:)'\newD1;
for i = 1:c2
    cogg1 =  cogg1 + ZCompFac(:,i)*sol_coeffs1(i); 
    %multiplying each element of the coefficient vector (alpha) to the
    %component and obtaining a linear combination for each sample
    cogg2 = cogg2 + ZCompFac(:,i)*sol_coeffs2(i);
end
CSC = [corr(ZPCs(:,1),cogg1)^2 corr(ZPCs(:,2),cogg2)^2]; %computing the squared correlation of COGG