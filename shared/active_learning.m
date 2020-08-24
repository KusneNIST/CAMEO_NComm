function [query, acc_ML, risks] = active_learning(nqueries, Y, W, initL)
% [query, acc_ML, risks] = active_learning(nqueries, Y, W, initL)
% Active learning on top of semi-supervised learning
% This is a new version which supports multiclass.
%
% Input:
%   nqueries: the desired number of active learning queries, excluding initL
%   Y: n*C matrix, each row in indicator encoding (one 1, others 0). The true labels of all points.
%      Note we don't require labeled points come first.
%   W: n*n graph weights matrix.  
%   initL: index between 1..n.  The initial labeled set.  
% Output:
%   The following are vectors of length nqueries.
%   query(i): new queries selected by active learning, given initL
%   acc_ML(i): The classification accuracy BEFORE adding query i
%   risks(i): expected risk of query(i), which is the smallest in all possible queries
%
%
% Please note there is no warranty.  In fact it is not even a software, but
% merely research code.  Feel free to modify it for your work.
% Xiaojin Zhu, zhuxj@cs.cmu.edu
% 2004


[n, nC] = size(Y);
L = initL;
U = setdiff(1:n, L);

% the laplacian
Delta = diag(sum(W))-W;

invDeltaU = full(inv(Delta(U,U))); % we don't expect it to be sparse anyway

% start active learning, find the best query to ask, add it to the training list, repeat.
for iteration=1:nqueries

	% compute f the harmonic solution.
	% f is a u*nC matrix
	f = invDeltaU * W(U,L) * Y(L,:);

	u = length(U);

	% The classification accuracy BEFORE adding new query
	[tmp ind] = sort(f, 2); predicted_class_U = ind(:,nC)-1;
	[tmp ind] = sort(Y(U,:), 2); true_class_U = ind(:,nC)-1;
	clear tmp
	acc_ML(iteration) = sum(true_class_U==predicted_class_U)/u;

	% nG(i,j) = G(i,j)/G(i,i), where G=invDeltaU
	%nG = invDeltaU ./ repmat(diag(invDeltaU)', u, 1);

	% we then compute f+(xk, yk), the harmonic function with one more 
	% labeled point (xk,yk).  

	% for efficiency, observe that 
	%     fplus = repmat(f(:,c), 1, u) + repmat(yk_c-f(:,c)', u, 1).*nG
	% can be decomposed into
	%     fplus = repmat(f(:,c), 1, u) - repmat(f(:,c)', u, 1).*nG + yk_c*nG
	% and we can precompute a maxfplus for all c=1:nC without the yk_c*nG term,
	% then for each yk=1:nC, we recompute fplus for c=yk and max it against the precomputed one
	% this is O(nC) instead of O(nC^2).
	for c=1:nC
		% fplus = repmat(f(:,c), 1, u) + repmat(-f(:,c)', u, 1).*nG;
		if (c==1)
			pre_maxfplus = repmat(f(:,c), 1, u) + repmat(-f(:,c)', u, 1).*(invDeltaU ./ repmat(diag(invDeltaU)', u, 1));
		else
			pre_maxfplus = max(pre_maxfplus, repmat(f(:,c), 1, u) + repmat(-f(:,c)', u, 1).*(invDeltaU ./ repmat(diag(invDeltaU)', u, 1)));
		end
	end
	risk = zeros(1, u); % the risk of querying point k
	for yk=1:nC
		c = yk; 
		% save some more space
		%fplus = repmat(f(:,c), 1, u) + repmat(1-f(:,c)', u, 1).*nG;
		%maxfplus = max(fplus, pre_maxfplus);
		%risk = risk + sum(1-maxfplus) .* f(:,yk)';
		risk = risk + sum(1- max(repmat(f(:,c), 1, u) + repmat(1-f(:,c)', u, 1).*(invDeltaU ./ repmat(diag(invDeltaU)', u, 1)) , pre_maxfplus)) .* f(:,yk)';
	end
	clear pre_maxfplus 

	% The following is old code with O(nC^2) time ------------------------
	%% This has to be repeated for yk = 1...nC
	%risk = zeros(1, u); % the risk of querying point k
	%for yk=1:nC
	%   
	%	% f is u*nC.  A 'slice' is a column vector for class c encoding.
	%	% For efficiency, we will compute f+(xk,yk)
	%	% one 'slice' (c) at a time, but for all candidate k:
	%	% f_plus is the 'slice', a u*u matrix, whose column k is [f+(xk,yk)]_class_c
	%	% 	f+(xk,yk)_c = f_c + (yk_c - fk_c)*nG
	%	% where yk_c is the c'th element for encoded yk.
	%	for c=1:nC
	%		if (yk==c) yk_c=1; else yk_c=0; end
	%		fplus = repmat(f(:,c), 1, u) + repmat(yk_c-f(:,c)', u, 1).*nG;

	%		% the risk for f+(xk,yk) is sum( 1 - max_c (f+(xk,yk)_c))
	%		if (c==1)
	%			maxfplus = fplus;
	%		else
	%			maxfplus = max(fplus, maxfplus);
	%		end
	%	end

	%	risk = risk + sum(1-maxfplus) .* f(:,yk)';
	%end
	% The above is old code with O(nC^2) time ------------------------

	% See which example in U, after being queried, will result in the minimum expected risk
	[minrisk, minUindex] = min(risk);

	query(iteration) = U(minUindex);
	risks(iteration) = minrisk;
	L = [L U(minUindex)];
	invDeltaU = compute_invAnegi(Delta(U,U), invDeltaU, minUindex);
	U = setdiff(U, U(minUindex));

end % iteration


