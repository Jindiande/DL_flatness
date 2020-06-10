clear all; clc;
% using TR method in paper
% "Complete Dictionary Recovery over the Sphere",
% to test faltness and error and there relation with noise level under 
% undercomplete dictionary case.

mu=1/100;
n=10% numberr of rows in dictionary
submatrix_num=10
n_sub=10% number of columns in dictioanry 
noise_level=[0:1e-1:5e-1]
sparsity=0.2;% about 20% of n
average_time=1
time_try=n_sub+20
p = 1e3;
tol = mu;

D_base=[]% n*n_sub baseline dictionary 
for i=1:(submatrix_num)
    D_base_sub=orth(randn(n_sub,n_sub))*diag([i:(n_sub+1-i)/(n_sub-1):n_sub+1]/(n_sub+1))*orth(randn(n_sub,n_sub))
    %D_base_sub=D_base_sub/norm(D_base_sub,'fro')% normalised Dictionary
    D_base = [D_base;D_base_sub]
end

%{
for ll = 1:p
   % generate p columns of sparse vectors for X
   X_base(randperm(n_sub,sparsity*n_sub),ll) = randn(sparsity*n_sub,1);
end
Y_base=D_base*X_base% baseline Y
%}

flatness_all_D=[]
error_all_X=[]
error_all_D=[]
size_all=[]
%{
X_noise_all=ones(length(noise_level),average_time,n_sub,p)
X_base_all=ones(length(noise_level),average_time,n_sub,p)
%}
for i=1:length(noise_level)
    
    %{
    
    Y_noise=Y_base+noise_level(i)*normrnd(0,1,[n,p])
    
    Y_noise_orth=real(inv(sqrtm(Y_noise*transpose(Y_noise)/(p*0.2)))*Y_noise)  %  orthogonal for Y_noise
    %Y_orth=inv(sqrtm(Y_base*transpose(Y_base)/(p*0.2)))*Y_base           % orthogonal for Y_base
    %}
    error_each_noise_X=[]
    error_each_noise_D=[]
    flat_each_noise_D=[]
    for k =1:(submatrix_num)
       error_sum_X=0
       error_sum_D=0
       flat_sum_D=0
       for time=1:average_time

            if_satis=1
            q_init=ones(n,1)
            D_noise_orth_inv_trans=[]% left inverse of dictionry
            [Y_noise,Y_base,X_base]=random_ini_Y(noise_level(i),D_base((k-1)*n_sub+1:k*n_sub,:),p,sparsity,n_sub)
            Y_noise_orth=real(inv(sqrtm(Y_noise*transpose(Y_noise)/(p*sparsity)))*Y_noise)  %  preconditioning for Y_noise
            for column=1:(time_try)% each time generate a column of left inverse of dictionary
                fprintf("i=%d,k=%d,time=%d,column=%d,real column%d\n",i,k,time,column,size(D_noise_orth_inv_trans,2))
                if size(D_noise_orth_inv_trans,2)~=0
                     q_init=null(D_noise_orth_inv_trans')*randn(n-size(D_noise_orth_inv_trans,2),1)% q_init lies in null space of D_noise
                     q_init=q_init/norm(q_init)
                else
                     q_init=randn(size(Y_noise_orth,1),1)
                     q_init=q_init/norm(q_init)
                end
                q = TR_Sphere(Y_noise_orth,mu,q_init,D_noise_orth_inv_trans); % using algorithm to generate q




                for j =1:size(D_noise_orth_inv_trans,2)% compare to other vector in res, if CosSim>0.95 ignore it.
                    if CosSim(q,D_noise_orth_inv_trans(:,j))>0.95 
                        if_satis=0
                        break
                    end

                end
                if if_satis==1
                   D_noise_orth_inv_trans=[D_noise_orth_inv_trans,q] 

                end
                if size(D_noise_orth_inv_trans,2)==n_sub-1
                   D_noise_orth_inv_trans=[D_noise_orth_inv_trans,null(transpose(D_noise_orth_inv_trans))] 
                   break
                end                
            end
            size_all=[size_all,size(D_noise_orth_inv_trans,2)]
            X_noise=transpose(D_noise_orth_inv_trans)*Y_noise_orth
            D_noise=sqrtm(Y_noise*transpose(Y_noise)/(p*sparsity))*transpose(D_noise_orth_inv_trans)% D_noise_orth_inv_trans is orthogonal
            sing_D=svd(D_noise)
            %D_noise=sqrtm(Y_noise*transpose(Y_noise)/(p*0.2))*inv(D_noise_orth_inv_trans*transpose(D_noise_orth_inv_trans))*D_noise_orth_inv_trans

            %sing_X=svd(X_noise)
            %flatness_sum=[flatness_sum;sing_X(1)/(sing_X(size(sing_X,1))+1e-10)]
            %X_noise_all(i,time,:,:)=X_noise
            %X_base_all(i,time,:,:)=X_base
            error_sum_X=error_sum_X+abs(norm(X_noise,'fro')-norm(X_base,'fro'))
            error_sum_D=error_sum_D+abs(norm(D_noise,'fro')-norm(D_base((k-1)*n_sub+1:k*n_sub,:),'fro'))
            flat_sum_D=flat_sum_D+sing_D(1)/(sing_D(size(sing_D,1))+1e-10)

       end
       error_each_noise_X=[error_each_noise_X;error_sum_X/average_time]
       error_each_noise_D=[error_each_noise_D;error_sum_D/average_time]
       flat_each_noise_D=[flat_each_noise_D,flat_sum_D/average_time]
    end

    %flatness_all=[flatness_all,flatness_sum]% flatness=min(singular value)/max(singular value)
    error_all_D=[error_all_D,error_each_noise_D]
    error_all_X=[error_all_X,error_each_noise_X]
    flatness_all_D=[flatness_all_D,flat_each_noise_D]
end

save('data3_complete')
%{
plot(noise_level,flatness_all-flatness_base)
xlabel('noise_level')
ylabel('flatness')
title('mean flatness of D in different noise level at')
legend('mean flateness')
figure(1)
%savefig(1,"D:\DL_code\simulation\TR_flatness_result\flatness_n_is"+int2str(n))
%}
plot(noise_level,error_all)
xlabel('noise_level')
ylabel('error')
title('mean error of X in different noise level at')
legend('error')
figure(1)
%savefig(1,"D:\DL_code\simulation\TR_flatness_result\error_n_is"+int2str(n))


