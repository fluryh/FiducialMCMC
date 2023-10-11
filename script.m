rng(1827)
num_iters = 5000;
burn_in=1000;
initial_guess = [.5,6];
n=100;
d=4;
data = zeros(d,n);
Sigma = buildSigma([.5,6],d);
for i =1:n
    data(:,i) = mvnrnd(zeros(d,1),Sigma);
end
[samples, accepts] = runConstrainedMH(num_iters,burn_in,@buildSigma,@DSigma, @isMA1param,initial_guess,data);
thetas = samples(1,:);
sigmas = samples(2,:);

max = -Inf;
MLEp = 0;
MLEs2 = 0;
for p = (1:100)/100
    for s2 = (1:1000)/100
        param = [p,s2];
        ll = ll_density(param,@buildSigma,data);
        if ll > max
            max = ll;
            MLEp = p;
            MLEs2 = s2;
        end
    end
end

histogram(sigmas)
xline(6, "r:","LineWidth",5.0)
xline(mean(sigmas),"g:","LineWidth",5.0)
xline(MLEs2, "b:", "LineWidth", 5.0)
legend("","True Value", "Simulated Mean", "MLE","FontSize",12)

histogram(thetas)
xline(.5, "r:","LineWidth",5.0)
xline(mean(thetas),"g:","LineWidth",5.0)
xline(MLEp, "b:", "LineWidth", 5.0)
legend("","True Value", "Simulated Mean", "MLE","FontSize",12)

disp(sum(accepts)/5000)