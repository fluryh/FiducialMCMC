function [samples, accepts] = runConstrainedMH(num_iters, burn_in, ...
    buildMatrixFn, SigmaDer, paramConstraints, ...
    initial_guess, data)
param_size = size(initial_guess,2);
d = size(data,1);
total_iters = num_iters + burn_in;
samples = zeros(param_size, total_iters+1);
accepts = zeros(1, total_iters+1);
samples(:,1) = initial_guess;
vec_size = size(initial_guess,1);
initial_Sigma = buildMatrixFn(initial_guess,d);
samp_num = 4; %this will need to be edited
%evens = [];
%odds = [];
%for i=1:(2^d)
    %bits = 2*int2bit(i,d)-1;
    %if prod(bits) == 1
        %evens = [evens, i];
    %else
        %odds = [odds,i];
    %end
%end
%Get the values for the initial guess so that it can compare
disp("Start SVD")
[S_init,D_init,~] = svd(initial_Sigma);
D_init = diag(sqrt(diag(D_init)));
disp("SVD Finished")
%if det(S_init) > 0
    %samp_array = evens(randperm(length(evens)));
%else
    %samp_array = odds(randperm(length(odds)));
%end
prev_jacProdSum = 0;
success_counter= 0;
counter = 1;
while success_counter < samp_num && counter <= 100
    bits = 2*randi(2,d-1,1) - 3;
    diagMat = diag([bits',(2*(det(S_init)>0)-1)*prod(bits)]);
    Stemp = S_init*diagMat;
    if rcond(eye(d) + Stemp) > 1e-10
        A = (eye(d) - Stemp)/(eye(d) + Stemp);
        A = (A - A')/2;
        %params_collapsed = collapseParameters(A,D_init);
        P = getProjectionMatrix(initial_guess,A,D_init,SigmaDer);
        prev_jacProdSum = prev_jacProdSum + getCoVJacobian(initial_guess,P,SigmaDer,A,D_init)*getFiducialJac(A,D_init,P,data,vec_size);
        %disp(getCoVJacobian(initial_guess,P,SigmaDer,A,D_init));
        %disp(getFiducialJac(A,D_init,P,data));
        success_counter = success_counter + 1;
     end
     disp(counter);
     counter = counter + 1;
end
prev_jacProdSum = prev_jacProdSum/success_counter;
prev_ll = ll_density(initial_guess,buildMatrixFn,data);
disp("Started!")
for iter=1:total_iters
    [proposal, flag] = generateProposal(samples(:,iter),paramConstraints);
    if flag == true
        disp("Couldn't find proposal")
        samples(:,iter+1) = samples(:,iter); 
    else
        %Get log-likelihood because this does not depend on signature
        %matrix
        log_likelihood = ll_density(proposal,buildMatrixFn,data);
        Sigma = buildMatrixFn(proposal,d);
        [S,D,~] = svd(Sigma);
        %if det(S) > 0
            %samp_array = evens(randperm(length(evens)));
        %else
            %samp_array = odds(randperm(length(odds)));
        %end
        jacProdSum = 0;
        success_counter= 0;
        counter = 1;
        while success_counter < samp_num %&& counter <= length(samp_array)
            bits = 2*randi(2,d-1,1) - 3;
            diagMat = diag([bits',(2*(det(S)>0)-1)*prod(bits)]);
            Stemp = S*diagMat;
            if rcond(eye(d) + Stemp) > 1e-10
                A = (eye(d) - Stemp)/(eye(d) + Stemp);
                A = (A - A')/2;
                %params_collapsed = collapseParameters(A,D);
                P = getProjectionMatrix(proposal,A,D,SigmaDer);
                jacProdSum = jacProdSum + getCoVJacobian(proposal,P,SigmaDer,A,D) * getFiducialJac(A,D,P,data,vec_size);
                %disp(getCoVJacobian(proposal,P,SigmaDer,A,D));
                %disp(getFiducialJac(A,D,P,data));
                success_counter = success_counter + 1;
            end
            counter = counter + 1;
        end
        jacProdSum = jacProdSum/success_counter;
        MH_ratio = log_likelihood + log(jacProdSum) - log(prev_jacProdSum) - prev_ll;
        %disp(MH_ratio);
        if log(rand) <= MH_ratio
            samples(:,iter+1)=proposal;
            accepts(:,iter)=1;
            prev_ll = log_likelihood;
            prev_jacProdSum=jacProdSum;
        else
            disp("Rejected by design")
            samples(:,iter+1)=samples(:,iter);
        end
    end
end
samples = samples(:,(burn_in+1):(num_iters+1));
accepts = accepts(:,(burn_in+1):(num_iters));
end