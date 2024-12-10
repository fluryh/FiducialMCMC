function Mat = dARIMA2(param,d)
phi = param(1);
theta = param(2);
sigma2 = param(3);
Mat = zeros(d*(d+1)/2+d,4);

ARMAvars = (1 + (phi+theta)^2/(1-phi^2));
ARMAvart = 2*sigma2*(phi+theta)/(1-phi^2);
ARMAvarp = 2*sigma2*((phi+theta)/(1-phi^2) + phi*(theta+phi)^2/((1-phi^2)^2));

ARMAcov = sigma2*(theta+phi+phi*((phi+theta)^2)/(1-phi^2));
ARMAcovs = (theta+phi+phi*(phi+theta)^2/(1-phi^2));
ARMAcovt = sigma2*(1+phi*(phi+2*theta)/(1-phi^2));
ARMAcovp = sigma2*(1+(phi+theta)^2/(1-phi^2) + ...
    phi*(2*phi+theta)/(1-phi^2) + 2*phi^2*(phi+theta)^2*((1-phi^2)^(-2)));
i = 1;
j= 1;
for index = 1:(d*(d+1)/2)
    for k = 1:i
        for l = 1:j
            if k==l
               Mat(index,1) = Mat(index,1) + ARMAvarp;
               Mat(index,2) = Mat(index,2) + ARMAvart;
               Mat(index,3) = Mat(index,3) + ARMAvars;
            else
               Mat(index,1) = Mat(index,1) + phi^(abs(k-l)-1) * ARMAcovp;
               Mat(index,1) = Mat(index,1) + (abs(k-l)-1)*phi^(abs(k-l)-2) * ARMAcov;
               Mat(index,2) = Mat(index,2) + phi^(abs(k-l)-1) * ARMAcovt;
               Mat(index,3) = Mat(index,3) + phi^(abs(k-l)-1) * ARMAcovs;
            end
        end
    end
    j = j+1;
    if j>d
        i = i+1;
        j=i;
    end
end
Mat(end-d+1:end,4) = 1;
end