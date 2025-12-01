function mat = DSigma(params,d)
theta = params(1);
sigma2 = params(2);
mat = zeros(d*(d+1)/2+d,2);
counter = 1;
for i=1:d
    for j=i:d
        if i==j
            mat(counter,:) = [1+theta^2,2*sigma2*theta];
        elseif j==i+1
            mat(counter,:) = [theta,sigma2];
        end
        counter = counter + 1;
    end
end
end