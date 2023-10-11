function P = getProjectionMatrix(proposal,A,D,dfn)
    d = size(D,1);
    derAL = zeros(d*(d+1)/2);
    I = eye(d);
    counter = 1;
    for i=1:d
        for j=(i+1):d
            J = zeros(d);
            J(j,i) = 1;
            J(i,j) = -1;
            B = (I+A) \ J / (I+A) *D^2 / (I-A) * (I+A);
            temp = 2*(B + B');
            counter2 = 1;
            for k = 1:d
                for l = k:d
                    derAL(counter,counter2) = temp(k,l);
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
                derAL(counter,counter2) = temp(k,l);
                counter2= counter2 + 1;
            end
        end
        counter = counter + 1;
    end
    derParam = dfn(proposal, d);
    M = derAL\derParam;
    P = M / (M' * M) * M';
    end