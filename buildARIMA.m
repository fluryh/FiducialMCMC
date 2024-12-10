function Mat = buildARIMA(param,d)
phi = param(1);
theta = param(2);
sigma2 = param(3);
Mat = diag(d);
var_help_array = zeros(1,d-1);
ARMAcov = sigma2*(theta+phi+phi*(phi+theta)^2/(1-phi^2));
Mat(1,1) = sigma2 * (1 + (phi+theta)^2/(1-phi^2));
for i = 1:d
    var_help_array(i) = phi^(i-1);
end
for i = 1:d
    if i ~= 1
        Mat(i,i) = i*Mat(1,1) + 2*ARMAcov*sum(var_help_array(1:(i-1)).*((i-1):-1:1));
    end
    for j = (i+1):d
        min_num = min(i,j-i);
        min_vec = min(min(1:(j-1),j-1:(j-1)),min_num);
        Mat(i,j) = i*Mat(1,1) + 2*ARMAcov*sum(var_help_array(1:(i-1)).*(i-(1:(i-1))));
        Mat(i,j) = Mat(i,j) + ARMAcov*sum(var_help_array(1:(j-1).*min_vec));
        Mat(j,i) = Mat(i,j);
    end
end
end