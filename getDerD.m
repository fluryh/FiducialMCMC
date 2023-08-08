function ders = getDerD(A,D,d)
    k = size(A,1);
    lambda = D(d,d);
    J = zeros(k);
    J(d,d) = 1;
    I = eye(k);
    ders = 2*lambda*(I-A)*inv(I+A)*J*inv(I-A)*(I+A);
end