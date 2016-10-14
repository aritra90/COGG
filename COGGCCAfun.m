function CCACOGGval = COGGCCAfun(ZPCs,ZCompFac)
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

[U_xx2,S_xx2,V_xx2] = svd(cov(ZPCs));
eig_XX2 = diag(S_xx2);
neweig_XX2 = 1./eig_XX2;
newSxx2 = diag(sqrt(neweig_XX2));
SigmaXX2 = U_xx2*newSxx2*V_xx2';

[U_yy,S_yy,V_yy] = svd(cov(ZCompFac)); 
eig_YY = diag(S_yy);
neweig_YY = 1./eig_YY; 
newSyy = diag(sqrt(neweig_YY));
SigmaYY = U_yy*newSyy*V_yy';

BigSigma2 = SigmaYY*covmat'*SigmaXX2;

[ls,ss,rs] = svd(BigSigma2);
d = diag(ss);
alpha = rs*SigmaXX2;
beta = SigmaYY*ls;
ccacogg1 = 0; 
ccacogg2 = 0;
coeffPCs = [];
coeffeverything = [];
for i = 1:length(rs)
    sum = alpha(i,1)*ZPCs(:,i);
    coeffPCs = [coeffPCs sum]; 
    ccacogg1 = ccacogg1 + sum;
end
for i = 1:length(ls)
    sum = beta(i,1)*ZCompFac(:,i);
    coeffeverything = [coeffeverything sum];
    ccacogg2 = ccacogg2 + sum;
end
CCACOGGval = corr(ccacogg1,ccacogg2)^2;
