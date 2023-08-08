function dc = DMA1LevelSet(param)
[A,D] = uncollapseParameters(param);
[n,~] = size(D);
dc = zeros(length(param)-2, length(param));
derCounter = 1;
for i=1:(n-2)
    for j=(i+2):n
        varCounter = 1;
        for k=1:(n-1)
            for l=(k+1):n
            derA = getDerA(A,D,k,l);
            dc(derCounter,varCounter) = derA(i,j);
            varCounter = varCounter + 1;
            end
        end
        for k=1:n
            derD = getDerD(A,D,k);
            dc(derCounter,varCounter) = derD(i,j);
            varCounter = varCounter + 1;
        end
        derCounter = derCounter + 1;
    end
end
for i = 1:(n-1)
    for j = i:(min(i+1,n-1))
        varCounter = 1;
        for k=1:(n-1)
            for l=(k+1):n
            derA = getDerA(A,D,k,l);
            dc(derCounter,varCounter) = derA(i,j) - derA(i+1,j+1);
            varCounter = varCounter + 1;
            end
        end
        for k=1:n
            derD = getDerD(A,D,k);
            dc(derCounter,varCounter) = derD(i,j) - derD(i+1,j+1);
            varCounter = varCounter + 1;
        end
        derCounter = derCounter + 1;
    end
end
end