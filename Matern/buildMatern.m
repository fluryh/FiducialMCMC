function Mat = buildMatern(param,d,dMat)
nu = param(1);
sigma2 = param(2);
rho = param(3);
Mat = diag(repelem(sigma2,d));
for i = 1:(d-1)
    for j = (i+1):d
        dist = dMat(i,j);
        Mat(i,j) = sigma2*2^(1-nu)/gamma(nu)*(sqrt(2*nu)*dist/rho)^nu*besselk(nu,sqrt(2*nu)*dist/rho);
        Mat(j,i) = Mat(i,j);
    end
end
end