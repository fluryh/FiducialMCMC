function Mat = buildMatern(param,d,dMat)
nu = param(1);
sigma2 = param(2);
rho = param(3);
Mat = diag(repelem(sigma2,d));
t1 = sigma2*2^(1-nu)/gamma(nu);
t2 = (sqrt(2*nu)/rho)^nu;
t3 = besselk(nu,sqrt(2*nu)*dist/rho);
for i = 1:(d-1)
    for j = (i+1):d
        Mat(i,j) = t1 * t2*(dMat(i,j)^nu)*t3;
        Mat(j,i) = Mat(i,j);
    end
end
end