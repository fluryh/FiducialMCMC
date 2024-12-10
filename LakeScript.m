rng(0191)
num_iters = 10000;
burn_in=0;
lake  = xlsread("LakeHuron.xlsx");
year = 1:98;
lm = fitlm(year,lake);
resid = lm.Residuals.Raw;
%resid = lake - mean(lake);
lakeARIMA = repelem(0,98);
lakeARIMA(1) = resid(1);
for i = 2:98
    lakeARIMA(i) = resid(i) + lakeARIMA(i-1);
end
initial_guess = [.25,.5,1,0];
step = [.01,.05,.1,.1];
tic;
[samples, accepts] = runConstrainedMH(num_iters,burn_in,@buildARIMA2,@meanARIMA,@dARIMA2,@isARIMAParam,initial_guess,step,lakeARIMA);
toc;
phis = samples(1,:);
thetas = samples(2,:);
sigmas = samples(3,:);
xnaughts = samples(4,:);