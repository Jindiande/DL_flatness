function [A,X, error,spar_pen] = learn_orthobasis_adm( Y_noise,Y, A_init, MAX_ITER, TOL, tau, DISPLAY) 

%%%
%
% learn_orthobasis_alternating
%
% Learns a sparsifying orthobasis
%
% Inputs:
%
% Outputs:
%
% Ju Sun and John Wright, January '14
%  Questions? johnwright@ee.columbia.edu
%  
% Last Update: Sat 18 Oct 2014 04:10:00 PM EDT 
%%%

done = false;
iter = 0;

dim = size(Y_noise, 1); 

A_old = A_init; 

while ~done,   
    
	iter = iter + 1;    
        
    X = prox_L1(A_old' * Y_noise, tau );    
    A_new = proj_orthogonal_group( Y_noise * X');         
    
    stepSize = norm(A_old - A_new, 'fro');
    
    if DISPLAY,
        disp(['Iteration ' num2str(iter) '  ||X||_1 ' num2str(norm1(X)) '  ||A-A_prev||_F ' num2str(stepSize)]);
    end

    if stepSize < TOL * sqrt(dim) || iter >= MAX_ITER,
    	done = true; 
    end 
    
    A_old = A_new; 
    
end

A = A_new; 

error = 0.5*norm(A'*Y - X,'fro')^2 ;

spar_pen=tau*norm(vec(X),1)
X=X
