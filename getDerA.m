function ders = getDerA(A,D,ARow,ACol)
    k = size(A,1);
    J = zeros(k);
    J(ACol,ARow) = 1;
    J(ARow,ACol) = -1;
    I = eye(k);
    B = I/(I+A)* J /(I+A)*D*D/(I-A)*(I+A);
    ders = 2 *(B + B.');
end