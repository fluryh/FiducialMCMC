function Sig = buildSigma(param, d)
sig2 = param(2);
theta = param(1);
Sig = zeros(d);
Sig(1,1) = sig2*(1+ theta^2);
Sig(1,2) = sig2*theta;
for i=2:(d-1)
    Sig(i,i-1) = sig2*theta;
    Sig(i,i) = sig2*(1+theta^2);
    Sig(i,i+1) = sig2*theta;
end
Sig(d,d) = sig2*(1+theta^2);
Sig(d,d-1) = sig2*theta;
end