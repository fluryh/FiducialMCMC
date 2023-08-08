function [A,D] = uncollapseParameters(param)
n = (sqrt(1 + 8 * size(param,2)) - 1)/2;
A = zeros(n);
counter = 1;
for i=1:(n-1)
    for j=(i+1):n
        A(i,j) = param(counter);
        A(j,i) = -1 * param(counter);
        counter = counter + 1;
    end
end
D = diag(param((end-n+1):end));
end