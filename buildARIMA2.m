function Mat = buildARIMA2(param,d)
phi = param(1);
theta = param(2);
sigma2 = param(3);
Mat = zeros(d);
ARMAcov = sigma2*(theta+phi+phi*((phi+theta)^2)/(1-phi^2));
ARMAvar = sigma2*(1 + (phi+theta)^2/(1-phi^2));
for i = 1:d
    for j = i:d
        for k = 1:i
            for l = 1:j
                if k==l
                    Mat(i,j) = Mat(i,j) + ARMAvar;
                else
                    Mat(i,j) = Mat(i,j) + phi^(abs(k-l)-1) * ARMAcov;
                end
            end
        end
        Mat(j,i) = Mat(i,j);
    end
end
end