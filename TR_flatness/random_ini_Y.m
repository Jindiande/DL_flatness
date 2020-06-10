function [Y_noise,Y_base,X_base]=random_ini_Y(noise_level,D_base,p,sparsity,n_sub)

for ll = 1:p
   % generate p columns of k-sparse vectors for Y
   X_base(randperm(n_sub,sparsity*n_sub),ll) = randn(sparsity*n_sub,1);
end


fprintf("%d %d\n",size(X_base,1),size(X_base,2))
Y_base=D_base*X_base
fprintf("%d %d\n",size(Y_base,1),size(Y_base,2))
Y_noise=Y_base+noise_level*normrnd(0,1,size(Y_base))% baseline Y