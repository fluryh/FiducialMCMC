function dARMAmat = dARMA(param,d)
phi = param(1);
theta = param(2);
sigma = param(3);
dARMAmat = zeros(d*(d+1)/2+d,4);
counter = 1;
for i = 1:d
    dARMAmat(counter,1) = 1 + ((theta + phi)^2)/(1-phi^2);
    dARMAmat(counter,2) = sigma * 2 * (theta + phi)/(1-phi^2);
    dARMAmat(counter,3) = sigma * 2 * (theta + phi)*(1+phi*theta)/(1-phi^2)^2;
    counter = counter + 1;
    for j = (i+1):d
        diff = max(0,j-i-1);
        dARMAmat(counter,1) = phi^diff * (theta + phi + phi * ((theta + phi)^2)/(1-phi^2));
        dARMAmat(counter,2) = phi^diff * sigma * (1 + 2*phi*(theta + phi)/(1-phi^2));
        if diff == 0
            dARMAmat(counter,3) = sigma * (phi^2*(theta^2+1)+ 4*phi*theta + theta^2 + 1)/(1-phi^2)^2;
        else
            dARMAmat(counter,3) = diff * phi^(diff-1) * sigma *(theta + phi + phi*((theta+phi)^2)/(1-phi^2)) ...
                + phi^diff * sigma * (phi^2*(theta^2+1)+4*phi*theta + theta^2 + 1)/(1-phi^2)^2;
        end
        counter = counter + 1;
    end
end
dARMAmat(end-d+1:end,4) = 1;
end