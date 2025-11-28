function derMat = dMaternMat(param,d,dMat)
nu = param(1);
sigma2 = param(2);
rho = param(3);
derMat = zeros(d*(d+1)/2+d,3);
counter = 1;
for i = 1:d
    for j = i:d
        if i == j
            derMat(counter,:) = [0,1,0];
        else
            dist = dMat(i,j);
            sigDer = 2^(1-nu)/gamma(nu)*(sqrt(2*nu)*dist/rho)^nu*besselk(nu,sqrt(2*nu)*dist/rho);
            nuDer = sigma2*2^(1-nu)/gamma(nu)*(sqrt(2*nu)*dist/rho)^nu*(besselk(nu,sqrt(2*nu)*dist/rho)*(-nu/2 - psi(nu) + log(sqrt(2*nu)*dist/rho)+1/(2*nu)) - .5*(sqrt(2*nu)*dist/rho)*(besselk(nu-1,sqrt(2*nu)*dist/rho)+besselk(nu+1,sqrt(2*nu)*dist/rho)));
            rhoDer = sigma2*2^(1-nu)/gamma(nu)*(sqrt(2*nu)*dist/rho)^nu*(-nu/rho*besselk(nu,sqrt(2*nu)*dist/rho) + sqrt(2*nu)*dist/(2*rho)*(besselk(nu-1,sqrt(2*nu)*dist/rho)+besselk(nu+1,sqrt(2*nu)*dist/rho)));
            derMat(counter,:) = [sigDer,nuDer,rhoDer];    
        end
        counter = counter + 1;
    end
end
end