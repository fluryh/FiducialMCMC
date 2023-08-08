function params = collapseParameters(A,D)
n = size(A,1);
params = 1:n*(n+1)/2;
counter = 1;
for i = 1:(n-1)
    for j = (i+1):n
        params(counter) = A(i,j);
        counter = counter + 1;
    end
end
for i =1:n
    params(counter) = D(i,i);
    counter = counter +1;
end
end