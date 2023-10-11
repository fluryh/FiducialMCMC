rng(1827)
num_iters = 5000;
burn_in=1000;
distmat = [1,1,1,sqrt(2); 1,1,sqrt(2),1; 1,sqrt(2),1,1;sqrt(2),1,1,1];
initial_guess = [.5,6,1];
n=100;
d=4;
data = zeros(d,n);
Sigma = buildMatern([.5,6,1],d,distmat);
for i =1:n
    data(:,i) = mvnrnd(zeros(d,1),Sigma);
end
buildFun = @(param,d)buildMatern(param,d,distmat);
dFun = @(param,d)dMaternMat(param,d,distmat);
[samples, accepts] = runConstrainedMH(num_iters,burn_in,buildFun,dFun,@isMaternParam,initial_guess,data);
nus = samples(1,:);
sigmas = samples(2,:);
rhos = samples(3,:);
