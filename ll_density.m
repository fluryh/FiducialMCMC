function logDens = ll_density(param,getCov,getMean,data)
[d,n] = size(data);
sigma = getCov(param,d);
mu = getMean(param,d);
logDens = -n * ((0.5*d)*log(2*pi)+(0.5)*log(det(sigma)));
for i=1:n
    logDens = logDens -0.5*((data(:,i)-mu')'/sigma*(data(:,i)-mu'));
end
end