function [U] = graph_cut_160402a(Y,U,S,w,B,X,C,gc_alpha)
% Y - spectra data for N samples, N x d
% U - membership matrix CxN
% S - sparse connectivity matrix, 'smoothness' matrix.
% w - graph cut weight multiplying the cost matrix.
% B - abundances
% X - endemembers
% C - composition
% gc_alpha - balance between two types of cost matricies.

% This software was developed by employees of the National Institute of
% Standards and Technology (NIST), an agency of the Federal Government and
% is being made available as a public service. Pursuant to title 17 United
% States Code Section 105, works of NIST employees are not subject to
% copyright protection in the United States.  This software may be subject
% to foreign copyright.  Permission in the United States and in foreign
% countries, to the extent that NIST may hold copyright, to use, copy,
% modify, create derivative works, and distribute this software and its
% documentation without fee is hereby granted on a non-exclusive basis,
% provided that this notice and disclaimer of warranty appears in all
% copies.

% THE SOFTWARE IS PROVIDED 'AS IS' WITHOUT ANY WARRANTY OF ANY KIND, EITHER
% EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, ANY
% WARRANTY THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED
% WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND
% FREEDOM FROM INFRINGEMENT, AND ANY WARRANTY THAT THE DOCUMENTATION WILL
% CONFORM TO THE SOFTWARE, OR ANY WARRANTY THAT THE SOFTWARE WILL BE ERROR
% FREE.  IN NO EVENT SHALL NIST BE LIABLE FOR ANY DAMAGES, INCLUDING, BUT
% NOT LIMITED TO, DIRECT, INDIRECT, SPECIAL OR CONSEQUENTIAL DAMAGES,
% ARISING OUT OF, RESULTING FROM, OR IN ANY WAY CONNECTED WITH THIS
% SOFTWARE, WHETHER OR NOT BASED UPON WARRANTY, CONTRACT, TORT, OR
% OTHERWISE, WHETHER OR NOT INJURY WAS SUSTAINED BY PERSONS OR PROPERTY OR
% OTHERWISE, AND WHETHER OR NOT LOSS WAS SUSTAINED FROM, OR AROSE OUT OF
% THE RESULTS OF, OR USE OF, THE SOFTWARE OR SERVICES PROVIDED HEREUNDER.

% A. Gilad Kusne, NIST, aaron.kusne@nist.gov, Release 8/01/2020
% If using this work for a publication, please cite:
% Kusne, A. Gilad, et al. "On-the-fly closed-loop materials discovery
% via Bayesian active learning." Nature communications 11.1 (2020): 1-11.

if(nargin > 4)
    if(gc_alpha < 1)
        U = graph_cut_advanced(Y,U,S,w,B,X,C,gc_alpha);
    end
else
    U = graph_cut_simple(Y,U,S,w);
end


% Simple Graph Cut --------------------------
function [U] = graph_cut_simple(Y,U,S,w)
[~,U]  = prune_empty_clusters_160402a(1,U);
N = size(Y,1);
k = size(U,1);
Dc = zeros([k N],'single');
[~, idx] = max(U);
c = zeros(k,size(Y,2));

% cluster mean
 for i=1:k
     if(sum(idx == i) > 1)
        c(i,:) = mean(Y(idx == i, :)); 
     else
        c(i,:) = Y(idx == i, :);
     end
 end

D = squareform(pdist([c; Y],'cosine'));
Dc = D(1:k,(k+1):end);
% if cost is NaN, set to maximum cost value.
Dc(isnan(Dc)) = max(Dc(:));

% smoothness cost
Sc = ones(k) - eye(k);

% run graphcut
gch = GraphCut('open', Dc*w, 10*Sc, S );
[gch L] = GraphCut('expand',gch);
gch = GraphCut('close', gch);

% convert cluster labels to cluster membership matrix
U = zeros(k,size(L,1));
for i=1:size(U,2)
    U(L(i)+1,i)=1;
end

% Advanced Graph Cut ------------------------------
function [U] = graph_cut_advanced(Y,U,S,w,B,X,C,gc_alpha)
N = size(Y,1);
k = size(U,1);

% initializat cost matrix
Dc = zeros([k N],'single');
[~, idx] = max(U);
c = zeros(k,size(Y,2));

% cluster mean
for i=1:k
    if (sum(idx == i) > 1)
        c(i,:) = mean(Y(idx == i, :)); 
    else
        c(i,:) = Y(idx == i,:);
    end
end

% cost Dc is computed as the distance of each spectra from cluster means.
D = squareform(pdist([c; Y],'cosine'));
Dc = D(1:k,(k+1):end);
Dc(isnan(Dc)) = max(Dc(:));

% cost Dist is computed as the sum residual of the difference between each
% spectra and its estimation using BX.
Dist = zeros(k,size(Y,1));
for i=1:k
    t = (Y - B{i}*X{i}).^2;
    Dist(i,:) = sum(t,2);
end
m = 0;
EPS = 0;

% normalize Dist by the sum of distances from all clusters.
SD      = sum(Dist,1);
Dc2      = (Dist)./repmat(SD, [k,1]);

% cost balance between the two types of costs
Dc = gc_alpha * Dc + (1-gc_alpha) * Dc2;

% smoothness cost
Sc = ones(k) - eye(k);

% compute graph cut
gch = GraphCut('open', Dc*w, 10*Sc, S );
[gch L] = GraphCut('expand',gch);
gch = GraphCut('close', gch);

% convert output to the cluster membership matrix
U = zeros(k,size(L,1));
for i=1:size(U,2)
    U(L(i)+1,i)=1;
end

'';
