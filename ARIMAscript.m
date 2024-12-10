rng(1827)
num_iters = 10000;
burn_in=1000;
n=100;
d=4;
initial_guess = [.4,.45,2,4];
step = 1.5*[.05,.05,.15,.15];
Sigma = buildARIMA2([.2,.25,4,3],d);
coverageProbs = [];
%mypool = parpool(10);
tic;
%for j= 1:200
    data = zeros(d,n);
    for i=1:n
        data(:,i) = mvnrnd(repelem(3,d),Sigma);
    end
    [samples, accepts] = runConstrainedMH(num_iters,burn_in,@buildARIMA2,@meanARIMA,@dARIMA2,@isARIMAParam,initial_guess,step,data);
    phis = samples(1,:);
    thetas = samples(2,:);
    sigmas = samples(3,:);
    xnaughts = samples(4,:);
    coverageProbs = [coverageProbs; [quantile(phis,.05),quantile(phis,.95),quantile(thetas,.05),quantile(thetas,.95),quantile(sigmas,.05),quantile(sigmas,.95),quantile(xnaughts,.05),quantile(xnaughts,.95)]];
%end
toc;
%delete(mypool);
%save("coverageProbSim")
%sum(coverageProbs(:,1) < .2 & coverageProbs(:,2) > .2)
