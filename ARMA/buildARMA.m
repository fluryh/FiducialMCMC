function ARMAmat = buildARMA(param,d)
phi = param(1);
theta = param(2);
sigma = param(3);
ARMAmat = sigma * (1 + ((theta + phi)^2)/(1-phi^2))*eye(d);
for dist = 1:(d-1)
    for j = 1:(d-dist)
        ARMAmat(j,j+dist) = sigma * phi^(dist-1) * (theta + phi + phi*((theta+phi)^2)/(1-phi^2));
        ARMAmat(j+dist,j) = ARMAmat(j,j+dist);
    end
end
end