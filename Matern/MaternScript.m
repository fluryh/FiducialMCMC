rng(1827)
addpath('..');
num_runs = 1;
MCMCEsts = zeros(3,num_runs);
MLEEsts = zeros(3,num_runs);
times = zeros(num_runs);
MCMCLB = zeros(3, num_runs);
MCMCUB = zeros(3, num_runs);
step = [.01,.1,.01];
num_iters = 20000;
burn_in=0;
distmat = ones(4,4);
distmat(1:size(distmat,1)+1:end) = 0;
distmat(1,2) = sqrt(5/4);
distmat(1,3) = sqrt(2);
distmat(2,3) = .5;
distmat(2,4) = sqrt(5/4);
distmat(tril(true(size(distmat)),-1)) = distmat(triu(true(size(distmat)),1));
n=100;
d=4;
buildFun = @(param,d)buildMatern(param,d,distmat);
dFun = @(param,d)dMaternMat(param,d,distmat);
likefun = @(param)-ll_density(param,buildFun,@meanMatern,data);
for run = 1:num_runs
    data = zeros(d,n);
    Sigma = buildFun([.5,6,1],d);
    for i =1:n
        data(:,i) = mvnrnd(zeros(d,1),Sigma);
    end
    tic;
    %profile on
    [samples, accepts] = runConstrainedMH(num_iters,burn_in,buildFun,...
        @meanMatern,dFun,@isMaternParam,initial_guess,step,data,8,.8);
    %profile viewer
    times(run) = toc;
    disp(times(run))
    disp(sum(accepts)/(num_iters-burn_in))
    MCMCEsts(:,run) = mean(samples,2);
    MCMCUB(:, run) = quantile(samples,.975,2);
    MCMCLB(:, run) = quantile(samples,.025,2);
    MLEEsts(:,run) = fminsearch(likefun,initial_guess);
end
%scatter(MCMCEsts(1,:), MCMCEsts(2,:))
%histogram(times)