rng(1828)
addpath('..');
num_runs = 1;
MCMCEsts = zeros(2,num_runs);
MLEEsts = zeros(2,num_runs);
times = zeros(num_runs);
MCMCLB = zeros(2, num_runs);
MCMCUB = zeros(2, num_runs);
step = [.01,.1];
num_iters = 100;
burn_in=0;
initial_guess = [.25,2];
n=4;
d=100;
likefun = @(param)-ll_density(param,@buildSigma,@meanMatern,data);
for run = 1:num_runs
    data = zeros(d,n);
    Sigma = buildSigma([.5,6],d);
    for i =1:n
        data(:,i) = mvnrnd(zeros(d,1),Sigma);
    end
    tic;
    %profile on
    [samples, accepts] = runConstrainedMH(num_iters,burn_in,@buildSigma, ...
        @meanMatern,@DSigma,@isMA1param,initial_guess,step,data,8,.8);
    %profile viewer
    times(run) = toc;
    disp(times(run))
    disp(sum(accepts)/(num_iters-burn_in))
    MCMCEsts(:,run) = mean(samples,2);
    MCMCUB(:, run) = quantile(samples,.975,2);
    MCMCLB(:, run) = quantile(samples,.025,2);
    MLEEsts(:,run) = fminsearch(likefun,[.25,2]);
end
%scatter(MCMCEsts(1,:), MCMCEsts(2,:))
%histogram(times)