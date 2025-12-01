function Mat = buildMatern(param,d,dMat)
nu = param(1);
sigma2 = param(2);
rho = param(3);
Mat = diag(repelem(sigma2,d));
t1 = sigma2*2^(1-nu)/gamma(nu) * (sqrt(2*nu)/rho)^nu;
t2 = sqrt(2*nu)/rho;
for i = 1:(d-1)
    for j = (i+1):d
        Mat(i,j) = t1 *(dMat(i,j)^nu)*besselk(nu,t2*dMat(i,j));
        Mat(j,i) = Mat(i,j);
    end
end
end