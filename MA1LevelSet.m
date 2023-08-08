function level_set_val = MA1LevelSet(param)
[A,D] = uncollapseParameters(param);
n = size(A,1);
tempZ = (eye(n) - A)/(eye(n)+A);
S = tempZ * D * D * tempZ'; 
level_set_val = [];
counter = 1;
for i = 1:(n-2)
    for j = (i+2):n
        level_set_val(counter) = S(i,j);
        counter = counter + 1;
    end
end
for i = 1:(n-1)
    for j = i:min(i+1,n-1)
        level_set_val(counter) = S(i,j) - S(i+1,j+1);
        counter = counter + 1;
    end
end
end