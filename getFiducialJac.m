function detJac = getFiducialJac(A,D,P,data,vec_size)
n = size(data,2);
d = size(data,1);
Jac = zeros(d*n,d*(d+1)/2);
[Q,~,~] = svds(P,vec_size);
counter = 1;
I = eye(d);
for i=1:d
    for j=(i+1):d
        J = zeros(d);
        J(i,j) = 1;
        J(j,i) = -1;
        temp = -2*I/(I + A) * J /(I-A) * data;
        for k=0:(n-1)
            Jac((d*k+1):(d*(k+1)),counter) = temp(:,k+1);
        end
        counter = counter + 1;
    end
end
for i=1:d
    J = zeros(d);
    J(i,i) = 1;
    temp = (I - A) /(I + A) * J /D *(I + A) /(I - A) * data;
    for k = 0:(n-1)
        Jac((d*k+1):(d*(k+1)),counter) = temp(:,k+1);
    end
    counter = counter + 1;
end
detJac = sqrt(det((Jac * Q)' * (Jac * Q)));
end