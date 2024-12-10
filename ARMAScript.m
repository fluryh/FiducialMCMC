rng(1827)
num_iters = 1000;
burn_in=0;
n=100;
d=4;
initial_guess = [.2,.25,.5,3];
step = [0.005,.005,.015,.015];
Sigma = buildARMA([.2,.25,.5,3],d);
data = zeros(d,n);
for i=1:n
    data(:,i) = mvnrnd(repelem(3,d),Sigma);
end
[samples, accepts] = runConstrainedMH(num_iters,burn_in,@buildARMA,@meanARMA,@dARMA,@isARMA11,initial_guess,step,data);
phis = samples(1,:);
thetas = samples(2,:);
sigmas = samples(3,:);
xnaughts = samples(4,:);
plot(phis)