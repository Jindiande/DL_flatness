% Using ADMM and preconditioning method to test relation between flatness of dictionary and noise of Y
% in Y=DX+E, D is undercomplete 
clear all; clc;

n=100% row number of Y
p=1e3% column number of Y
n_sub=10% row number of submatrix
sparsity=0.2% sparsity of X_base
noise_level=[0:5e-2:5e-1]

MAX_ITER = 10000 % 
average_time=5
DISPLAY  = false 
TOL = 1e-5;
%tau = [0.1:0.1:2];
tau=0.1% fix tau to 0.1
% baseline A where flatness of each submatrix range from 1/n_sub+1 to n_sub/n_sub+1
flatness_base=[1/(n_sub+1):1/(n_sub+1):n_sub/(n_sub+1)]
D_base=[]% n*n_sub baseline dictionary 

for i=1:(n/n_sub)
    D_base_sub=orth(randn(n_sub,n_sub))*diag([i:(n_sub+1-i)/(n_sub-1):n_sub+1])*orth(randn(n_sub,n_sub))
    D_base_sub=D_base_sub/norm(D_base_sub,'fro')% normalised Dictionary
    D_base = [D_base;D_base_sub]
end

%{
for ll = 1:p
   % generate p columns of k-sparse vectors for Y
   X_base(randperm(n_sub,sparsity*n_sub),ll) = randn(sparsity*n_sub,1);
end
Y_base=D_base*X_base% baseline Y
%}
flatness_all=[]
error_all=[]

for i =1:length(noise_level)
    %{
    Y_noise=Y_base+noise_level(i)*normrnd(0,1,[n,p])
    Y_noise_orth=inv(sqrtm(Y_noise*transpose(Y_noise)/(p*0.2)))*Y_noise  %  orthogonal for Y_noise
    Y_orth=inv(sqrtm(Y_base*transpose(Y_base)/(p*0.2)))*Y_base           % orthogonal for Y_base
    %}
    flatness_each_noise=[]
    error_each_noise=[]
    for j=1:(n/n_sub)
        flatness_sum=0
        error_sum=0
        flatness_ave=0
        error_ave=0
        for time =1:average_time
            [Y_noise,Y_base]=random_ini_Y(noise_level(i),D_base,n,p,sparsity,n_sub)
            Y_noise_orth=real(inv(sqrtm(Y_noise*transpose(Y_noise)/(p*sparsity))))*Y_noise  %  orthogonal for Y_noise
            Y_orth=real(inv(sqrtm(Y_base*transpose(Y_base)/(p*sparsity))))*Y_base           % orthogonal for Y_base
           [D_orth_sub,X,error,spar_pen]=learn_orthobasis_adm( Y_noise_orth((j-1)*n_sub+1:j*n_sub,:),Y_orth((j-1)*n_sub+1:j*n_sub,:), proj_orthogonal_group(randn(n_sub,n_sub)), MAX_ITER, TOL, tau, DISPLAY);
            D_sub=sqrtm(Y_noise((j-1)*n_sub+1:j*n_sub,:)*transpose(Y_noise((j-1)*n_sub+1:j*n_sub,:))/(p*sparsity))*D_orth_sub
            sing_D=svd(D_sub)
            %sing_D_base=svd(D_base((j-1)*n_sub+1:j*n_sub,:))
            error_sum=error_sum+abs(norm(D_sub,'fro')-norm(D_base((j-1)*n_sub+1:j*n_sub,:),'fro'))% distance between eigenvalues of D_base and D_result
            flatness_sum=flatness_sum+sing_D(n_sub)/(sing_D(1)+1e-10)
        end
        flatness_ave=flatness_sum/average_time
        error_ave=error_sum/average_time
        flatness_each_noise=[flatness_each_noise;flatness_ave]   
        error_each_noise=[error_each_noise;error_ave]
    end
    flatness_all=[flatness_all,flatness_each_noise]%10*i
    error_all=[error_all,error_each_noise]%10*i

end
save('june10')
%{
plot(noise_level,mean(flatness_all,1))
xlabel('noise_level')
ylabel('flatness')
title('mean flatness in different noise level at')
legend('mean flateness')
figure(1)
savefig(1,"ADM_flatness_result\flatness_noise_Dic\flatness_n_is"+int2str(n))
%}

Legend=cell(10,1)
for i = 1:10
    plot(noise_level,flatness_all(i,:)-flatness_base(i),'-o')
    Legend{i}=strcat('No', num2str(i),'sub_matrix');
    hold on
end
hold off
legend(Legend)
title('flatness in different noise level at')
figure(1)
savefig(1,"D:\DL_code\image_experiment_focm\ADM_flatness_result\flatness_noise_Dic\flatness_n_is"+int2str(n))
%{
plot(noise_level,mean(error_all,1))
xlabel('noise_level')
ylabel('error')
title('mean flatness in different noise level at')
legend('mean flateness')
figure(1)
savefig(1,"ADM_flatness_result\flatness_noise_Dic\error_n_is"+int2str(n))
%}

Legend=cell(10,1)
for i = 1:10
    plot(noise_level,error_all(i,:),'-o')
    Legend{i}=strcat('No', num2str(i),'sub_matrix');
    hold on
end
hold off
legend(Legend)
title('error in different noise level at')
figure(1)
savefig(1,"D:\DL_code\image_experiment_focm\ADM_flatness_result\flateness_noise_Dic\error_n_is"+int2str(n))

function [Y_noise,Y_base]=random_ini_Y(noise_level,D_base,n,p,sparsity,n_sub)

for ll = 1:p
   % generate p columns of k-sparse vectors for Y
   X_base(randperm(n_sub,sparsity*n_sub),ll) = randn(sparsity*n_sub,1);
end


fprintf("%d %d\n",size(X_base,1),size(X_base,2))
Y_base=D_base*X_base
Y_noise=Y_base+noise_level*normrnd(0,1,[n,p])% baseline Y
end