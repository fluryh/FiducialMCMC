function [samples, accepts] = runConstrainedMH(num_iters, burn_in, ...
    buildMatrixFn, buildMeanFn, SigmaDer, paramConstraints, ...
    initial_guess, step_size, data)
param_size = size(initial_guess,2);
d = size(data,1);
n = size(data,2);
total_iters = num_iters + burn_in;
samples = zeros(param_size, total_iters+1);
accepts = zeros(1, total_iters+1);
samples(:,1) = initial_guess;
initial_Sigma = buildMatrixFn(initial_guess,d);
initial_mu = buildMeanFn(initial_guess,d);
samp_num = 8; %This will need to be edited
keep_prop = .8; %This will need to be edited
nextAs = zeros(samp_num,d);
[S_init,D_init,~] = svd(initial_Sigma);
D_init = diag(sqrt(diag(D_init)));
prev_jacProdSum = 0;
success_counter= 0;
%mypool = parpool(8);
for iter=1:samp_num
    bits = 2*randi(2,d-1,1) - 3;
    bits = [bits',(2*(det(S_init)>0)-1)*prod(bits)];
    diagMat = diag(bits);
    Stemp = S_init*diagMat;
    nextAs(success_counter+1,:) = bits;
    if rcond(eye(d) + Stemp) > 1e-5
        A = (eye(d) - Stemp)/(eye(d) + Stemp);
        A = (A - A')/2;
        Jac = zeros(d*n,d*(d+1)/2+d);
        counter = 1;
        I = eye(d);
        for i=1:d
            for j=(i+1):d
                J = zeros(d);
                J(i,j) = 1;
                J(j,i) = -1;
                temp = -2*I/(I + A) * J /(I-A) * (data-initial_mu');
                for k=0:(n-1)
                    Jac((d*k+1):(d*(k+1)),counter) = temp(:,k+1);
                end
                counter = counter + 1;
            end
        end
        for i=1:d
            J = zeros(d);
            J(i,i) = 1;
            temp = (I - A) /(I + A) * J /D_init *(I + A) /(I - A) * (data-initial_mu');
            for k = 0:(n-1)
                Jac((d*k+1):(d*(k+1)),counter) = temp(:,k+1);
            end
        counter = counter + 1;
        end
        for k = 0:(n-1)
            Jac((d*k+1):(d*(k+1)),((d*(d+1)/2)+1):(d*(d+1)/2 + d)) = eye(d);
        end
        derAL = zeros(d*(d+1)/2+d);
        counter = 1;
        for i=1:d
            for j=(i+1):d
                J = zeros(d);
                J(j,i) = 1;
                J(i,j) = -1;
                B = (I+A) \ J / (I+A) *D_init^2 / (I-A) * (I+A);
                temp = 2*(B + B');
                counter2 = 1;
                for k = 1:d
                    for l = k:d
                        derAL(counter,counter2) = temp(k,l);
                        counter2 = counter2 + 1;
                    end
                end
                counter = counter + 1;
            end
        end
        for i =1:d
            J = zeros(d);
            J(i,i) = 1;
            temp = 2 * D_init(i,i) * (I-A) / (I+A) * J / (I-A) * (I+A);
            counter2 = 1;
            for k =1:d
                for l = k:d
                    derAL(counter,counter2) = temp(k,l);
                    counter2= counter2 + 1;
                end
            end
            counter = counter + 1;
        end
        derAL(end-d+1:end,end-d+1:end) = eye(d);
        derParam = SigmaDer(initial_guess, d);
        JacMat = Jac/derAL*derParam;
        prev_jacProdSum = prev_jacProdSum + sqrt(det(JacMat' * JacMat));
        success_counter = success_counter + 1;
    end
end
prev_jacProdSum = prev_jacProdSum/success_counter;
prev_ll = ll_density(initial_guess,buildMatrixFn,buildMeanFn,data);
for iter=1:total_iters
    if iter == 313
        1 + 1;
    end
    [proposal, flag] = generateProposal(samples(:,iter),paramConstraints,step_size);
    if flag == true
        %disp("Couldn't find proposal")
        samples(:,iter+1) = samples(:,iter); 
    else
        log_likelihood = ll_density(proposal,buildMatrixFn,buildMeanFn,data);
        Sigma = buildMatrixFn(proposal,d);
        mu = buildMeanFn(proposal,d);
        [S,D,~] = svd(Sigma);
        derParam = SigmaDer(proposal, d);
        jacProdSum = 0;
        success_counter= 0;
        keepAs = nextAs;
        nextAs = keepAs(randperm(size(keepAs,1)),:);
        for iter_counter=1:samp_num
            if  iter_counter <= keep_prop*samp_num
                bits = nextAs(iter_counter,:);
            else
                bits = 2*randi(2,d-1,1) - 3;
                bits = [bits',(2*(det(S)>0)-1)*prod(bits)];
                nextAs(iter_counter,:) = bits;
            end
            diagMat = diag(bits);
            Stemp = S*diagMat;
            if rcond(eye(d) + Stemp) > 1e-5
                I = eye(d);
                A = (I - Stemp)/(I + Stemp);
                A = (A - A')/2;
                Jac = zeros(d*n,d*(d+1)/2+d);
                counter = 1;
                tempMat1 = -2*I/(I + A);
                tempMat2 = (I-A) \ (data-mu');
                for i=1:d
                    for j=(i+1):d
                        J = zeros(d);
                        J(:,j) = -tempMat1(:,i);
                        J(:,i) = tempMat1(:,j);
                        temp = J * tempMat2 ;
                        for k=0:(n-1)
                            Jac((d*k+1):(d*(k+1)),counter) = temp(:,k+1);
                        end
                        counter = counter + 1;
                    end
                end
                tempMat1 = (I - A) /(I + A);
                tempMat2 = D \(I + A) /(I - A) * (data-mu');
                for i=1:d
                    J = zeros(d);
                    J(:,i) = tempMat1(:,i);
                    temp = J * tempMat2;
                    for k = 0:(n-1)
                        Jac((d*k+1):(d*(k+1)),counter) = temp(:,k+1);
                    end
                counter = counter + 1;
                end
                for k = 0:(n-1)
                    Jac((d*k+1):(d*(k+1)),((d*(d+1)/2)+1):(d*(d+1)/2 + d)) = eye(d);
                end
                derAL = zeros(d*(d+1)/2+d);
                counter = 1;
                BMat1 = I/(I+A);
                BMat2 = (I+A)\ D^2 / (I-A) * (I+A);
                for i=1:d
                    for j=(i+1):d
                        J = zeros(d);
                        J(:,j) = -BMat1(:,i);
                        J(:,i) = BMat1(:,j);
                        B1 =  J * BMat2;
                        temp1 = 2*(B1 + B1');
                        counter2 = 1;
                        for k = 1:d
                            for l = k:d
                                derAL(counter,counter2) = temp1(k,l);
                                counter2 = counter2 + 1;
                            end
                        end
                        counter = counter + 1;
                    end
                end
                tempMat1 = (I-A) / (I+A);
                tempMat2 = (I-A) \ (I+A);
                for i =1:d
                    J = zeros(d);
                    J(:,i) = tempMat1(:,i);
                    temp2 = 2 * D(i,i) * J * tempMat2 ;
                    counter2 = 1;
                    for k =1:d
                        for l = k:d
                            derAL(counter,counter2) = temp2(k,l);
                            counter2= counter2 + 1;
                        end
                    end
                    counter = counter + 1;
                end
                derAL(end-d+1:end,end-d+1:end) = eye(d);
                JacMat = Jac/derAL * derParam;
                %disp(sqrt(det(JacMat' * JacMat)))
                jacProdSum = jacProdSum + sqrt(det(JacMat' * JacMat));
                success_counter = success_counter + 1;
            end
        end
        jacProdSum = jacProdSum/success_counter;
        MH_ratio = log_likelihood + log(jacProdSum) - log(prev_jacProdSum) - prev_ll;
        disp(log_likelihood)
        disp(prev_ll)
        disp(log(jacProdSum))
        disp(log(prev_jacProdSum))
        if log(rand) <= MH_ratio
            samples(:,iter+1)=proposal;
            accepts(:,iter)=1;
            prev_ll = log_likelihood;
            prev_jacProdSum=jacProdSum;
        else
            samples(:,iter+1)=samples(:,iter);
            nextAs = keepAs;
        end
    end
end
samples = samples(:,(burn_in+1):(num_iters+1));
accepts = accepts(:,(burn_in+1):(num_iters));
%delete(mypool);
end