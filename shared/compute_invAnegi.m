function invAnegi = compute_invAnegi(A, invA, i)
% invAnegi = compute_invAnegi(A, invA, i)
%
% Compute inv(A_neg_i), where A_neg_i is the matrix by removing the ith row/column of A.
%
% A: a n*n matrix.
% invA: A^{-1}, a n*n matrix.
% i: the row/column to be removed.
%
% invAnegi: the (n-1)*(n-1) matrix inv(A_neg_i).
%
% Algorithm: see section 5 of my active learning note.

n = size(invA,1);

% create a perm = (i,1,2,3,...,i-1, i+1, ...)
perm=[i 1:i-1 i+1:n]; 

B = A(perm,perm);
invB = invA(perm, perm);

u=[-1 zeros(1,n-1)]';
v=(B(1,:) - [1 zeros(1,n-1)])';

tmp = v' * invB;
invBprime = invB + invB(:,1) * tmp / (1 - v' * invB(:,1));
%invBprime = invB - invB * u * v' * invB / (1 + v' * invB * u);

w = B(:,1); w(1)=0;
tmp = invBprime * w;
invBprimeprime = invBprime + tmp * invBprime(1,:) / (1 - invBprime(1,:) * w);
%invBprimeprime = invBprime - invBprime * w * u' * invBprime / (1 + u' * invBprime * w);

invAnegi = invBprimeprime(2:n,2:n);
