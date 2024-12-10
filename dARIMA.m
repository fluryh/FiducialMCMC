function Mat = dARIMA(param,d)
phi = param(1);
theta = param(2);
sigma2 = param(3);
[~,param_dim] = size(param);
Mat = zeros(d*(d+1)/2,param_dim);
var_help_array = zeros(1,d-1);
dvar_help_array = zeros(1,d-1);
ARMAcov = sigma2*(theta+phi+phi*(phi+theta)^2/(1-phi^2));
ARMAcovs = (theta+phi+phi*(phi+theta)^2/(1-phi^2));
ARMAcovt = sigma2*(1+phi*(phi+2*theta)/(1-phi^2));
ARMAcovp = sigma2*(1+(phi+theta)^2/(1-phi^2) + ...
    phi*(2*phi+theta)/(1-phi^2) + 2*phi^2*(phi+theta)^2*((1-phi^2)^(-2)));
Mat(1,1) = sigma2 * ((2*phi+theta)/(1-phi^2) + 2*phi*(phi+theta)^2/((1-phi^2)^2));
Mat(1,2) = sigma2 * (2*theta+phi)/(1-phi^2);
Mat(1,3) = (1 + (phi+theta)^2/(1-phi^2));
for i = 1:d
    var_help_array(i) = phi^(i-1);
    dvar_help_array(i) = (i-1)*phi^(i-2);
end
j = 2;
i = 1;
for index = 2:(d*(d+1)/2)
    if i==j
        Mat(index,1) = i*Mat(1,1) + 2*ARMAcov*sum(dvar_help_array(1:(i-1)).*((i-1):-1:1))...
            + 2*ARMAcovp*sum(var_help_array(1:(i-1)).*((i-1):-1:1));
        Mat(index,2) = i*Mat(1,2) + 2*ARMAcovt*sum(var_help_array(1:(i-1)).*((i-1):-1:1));
        Mat(index,3) = i*Mat(1,3) + 2*ARMAcovs*sum(var_help_array(1:(i-1)).*((i-1):-1:1));
    else
        min_num = min(i,j-i);
        min_vec = min(min(1:(j-1),j-1:(j-1)),min_num);
        Mat(index,1) = i*Mat(1,1) + 2*ARMAcovp*sum(var_help_array(1:(i-1)).*(i-(1:(i-1))))...
            + 2*ARMAcov*sum(dvar_help_array(1:(i-1)).*(i-(1:(i-1)))) + ARMAcovp*sum(var_help_array(1:(j-1).*min_vec))...
            + ARMAcov*sum(dvar_help_array(1:(j-1).*min_vec));
        Mat(index,2) = i*Mat(1,2) + 2*ARMAcovt*sum(var_help_array(1:(i-1)).*(i-(1:(i-1)))) + ARMAcovt*sum(var_help_array(1:(j-1).*min_vec));
        Mat(index,3) = i*Mat(1,3) + 2*ARMAcovs*sum(var_help_array(1:(i-1)).*(i-(1:(i-1)))) + ARMAcovs*sum(var_help_array(1:(j-1).*min_vec));
    end
    j = j+1;
    if j > d
        i = i+1;
        j = i;
    end
end
end