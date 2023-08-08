function CoVdetJac = getCoVJacobian(params,P,paramDer,A,D)
d = size(A,1);
[Q,~,~] = svds(P,size(params,2));
detMat = zeros(d*(d+1)/2);
I = eye(d);
counter = 1;
for i=1:d
    for j=(i+1):d
        J = zeros(d);
        J(j,i) = 1;
        J(i,j) = -1;
        B = inv(I + A) * J / (I+A) *D^2 / (I-A) *(I +A);
        temp = 2*(B + B');
        counter2 = 1;
        for k = 1:d
            for l = k:d
                detMat(counter,counter2) = temp(k,l);
                counter2 = counter2 + 1;
            end
        end
        counter = counter + 1;
    end
end
for i =1:d
    J = zeros(d);
    J(i,i) = 1;
    temp = 2 * D(i,i) * (I-A) / (I+A) * J / (I-A) * (I+A);
    counter2 = 1;
    for k =1:d
        for l = k:d
            detMat(counter,counter2) = temp(k,1);
            counter2= counter2 + 1;
        end
    end
    counter = counter + 1;
end
det2 = det((Q'* detMat') * (detMat * Q));
derivative = paramDer(params,d);
det1 = det(derivative'*derivative);
CoVdetJac = det1/det2;
end