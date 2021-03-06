function [x,optval] = Solve_TR_Subproblem(A,b,Delta)

% solve_trust_region_subproblem by SDP relaxation
%   min (1/2) z' A z + z' b   s.t.  ||z||_2 <= Delta
%
% for z \in R^n via an (exact) convex relaxation by SDP lifting
%   min   < [ A b; b' 0 ], Z >  s.t. Z >= 0, Z_{n+1,n+1} = 1, trace(Z) = 2.

	G = (1/2) * [ A, b; b', 0 ];
	n = length(b);

	cvx_quiet(1);

	cvx_begin sdp
		variable X(n+1,n+1) symmetric;
		minimize(trace(G'*X))
		subject to
		trace(X) <= 1+Delta^2;
		X(n+1,n+1) == 1;
		X >= 0;
	cvx_end

	% recover the optimal direction by the leading eigenvector using SVD 
	[u,~,~] = svds(X,1);
	x = u(1:n);

	optval = cvx_optval;
	nx = - x;
	if (1/2) * nx'*A*nx + b'*nx < (1/2) * x'* A * x + b' * x
		x = nx;
	end

end
