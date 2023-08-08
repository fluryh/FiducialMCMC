function logDens = ll_density(param,getCov,data)
[d,n] = size(data);
sigma = getCov(param,d);
logDens = 0;
for i=1:n
    logDens = logDens - (0.5*d)*log(2*pi)-(0.5)*log(det(sigma))-0.5*(data(:,i)'/sigma*data(:,i));
end
end