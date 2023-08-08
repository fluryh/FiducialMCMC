num_iters = 5000;
burn_in=100;
initial_guess = [.75,2];
n=1;
d=100;
data = zeros(d,n);
Sigma = buildSigma([.5,6],d);
for i =1:n
    data(:,i) = mvnrnd(zeros(d,1),Sigma);
end
[samples, accepts] = runConstrainedMH(num_iters,burn_in,@buildSigma,@DMA1LevelSet, @DSigma, @isMA1param,initial_guess,data);
thetas = samples(1,:);
sigmas = samples(2,:);
plot(sigmas)
plot(thetas)