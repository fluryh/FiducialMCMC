function valid = isMA1param(params)
%params assumed to be rho then sigma
valid = abs(params(1)) <= 1 && params(2) > 0; 
end