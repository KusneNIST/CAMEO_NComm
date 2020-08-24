function [E_,P_,U_,Cond_,numclust, iter_breaks]=Phase_Mapping_160715a(data,parameters)

% This function computes multi-model endmembers and their respective abundances
% Input:
%   - X         : Data (N x D) matrix. N data points of dimensionality D. 
%   - parameters: The parameters set by PCOMMEND_Parameters function.
%
% Output:
%   - E         : Cell of C endmembers matrices. One MxD matrix per cluster.
%   - P         : Cell of C abundance matrices. One NxM matrix per cluster.
%   - U         : Fuzzy membership matrix CxN.
%
%
% ORIGINAL CODE COPYRIGHT:
% This product is Copyright (c) 2013 University of Missouri, University
% of Florida and University of Louisville.
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions
% are met:
%
%   1. Redistributions of source code must retain the above copyright
%      notice, this list of conditions and the following disclaimer.
%   2. Redistributions in binary form must reproduce the above copyright
%      notice, this list of conditions and the following disclaimer in the
%      documentation and/or other materials provided with the distribution.
%   3. Neither the name of the University nor the names of its contributors
%      may be used to endorse or promote products derived from this software
%      without specific prior written permission.
%
% MODIFICATIONS COPYRIGHT:
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
% Kusne, A. G., et al. "On-the-fly Closed-loop Autonomous Materials
% Discovery via Bayesian Active Learning." arXiv preprint arXiv:2006.06141
% (2020).
%%

% Initialize parameters
% me: removed U from output
% C = Number of clusters and number of endmember sets.
% M = Number of Endmembers per endmember set. 
% X = Cell of C endmembers matrices. One (Mxd) matrix per cluster.
% X{i} initialized to matrix with first M spectra in cluster i.
% U = initialized by input;
% PAST: U = the fuzzy membership matrix (CxN), set by PCOMMEND_fcm(X,C) = Fuzzy C Means
% PAST: http://en.wikipedia.org/wiki/Fuzzy_clustering#Fuzzy_c-means_clustering
% [E,U] = PCOMMEND_Initialize(X,parameters.C,parameters.M);

% -------- Unpack Data -------------------------------
% (X,parameters,U,S,xy,eCMP,TTH)
X = data.spectra;
U = data.spectra_cluster_labels;
S = data.connectivity;
xy = data.coordinates;
eCMP = data.composition;
TTH = data.spectra_sampling_coordinates;
w = parameters.graphcut.wgc;
gc_alpha = parameters.graphcut.alpha;
parameters = parameters.PCOMMEND;
[~,idx1] = max(U);
E_ = {};
P_ = {};
U_ = {};
C_ = [];
Cond_ = [];
iter_breaks = {};
numclust = 0;

%--------Initialization -------------------------
U = graph_cut_160402a(X,U,S,w); % graph cut the data
[parameters.C,U]  = prune_empty_clusters_160402a(1,U); % remove empty clusters
if(sum((sum(U,2)>0))<=1)
 %   'only 1'
    numclust = 1;    
else
% [~,idx2] = max(U);   
E = initialize_E(X,U,parameters.M);
Num_Cluster_Going_In = parameters.C; %length(unique(idx2));

[N,D] = size(X); 
Cond_old = inf;


sc = sum(U,2);
for i = 1:Num_Cluster_Going_In
    tm = min([parameters.M sc(i)]);
      P{i} = ones(N, tm)*(1/tm) ;
end
Num_Cluster_Current = Num_Cluster_Going_In;

for iter = 1:parameters.iterationCap 
    
    %Update abundances B ------------------------
    P_old = P;
    P     = PCOMMEND_P_update(X,E,Num_Cluster_Current,parameters.EPS);
        
    %Update Endmembers E ------------------------
    E_old = E;
    % m = Fuzzifier for membership.
    % alpha = Regularization Parameter. Range of alpha should be between 0 and 1.
    % Equation 6 for updating endmembers based on objective function.
    E      = PCOMMEND_E_Update_NonNegative(U,parameters.m,parameters.alpha,P,X,parameters.EPS);
        
    % Update membership U ----------------------
    U_old = U;
    U = graph_cut_160402a(X,U,S,w,P,E,eCMP,gc_alpha);
    
    min_cluster_size = 1;
	[parameters.C,U,E,P,U_old,E_old,P_old]  = prune_empty_clusters_160402a(min_cluster_size,U,E,P,U_old,E_old,P_old);
        
    % Kill this iteration if the number of clusters has changed.
    [~,idx4] = max(U);
    Num_Cluster_Current = length(unique(idx4));
    
    % breaking conditions for ending a rapid drop in number of clusters.
    if iter > 1
        if(Num_Cluster_Current < parameters.GRENDEL_lower_bound_cutoff * C_(1))
            P_ = P_(1);
            E_ = E_(1);
            U_ = U_(1); 
            numclust = numclust(1);
            Cond_ = Cond_(1);
            iter_breaks{1} = [iter Num_Cluster_Going_In Num_Cluster_Current];
            break
        end
    end
    
    P_{iter} = P;
    E_{iter} = E;
    U_{iter} = U;
    C_ = [C_ Num_Cluster_Current];
    
    % Stopping criteria ------------------------------
    % minimizing change in U, P, and E. (membership, proportion, and endmembers).
    Cond  = norm(U - U_old);
    for i = 1:Num_Cluster_Current
        Cond = Cond + norm(P{i} - P_old{i}) + norm(E{i} - E_old{i}) ;
    end
    Cond_(iter) = Cond;
    if(abs(Cond - Cond_old) <parameters.changeThresh)   
        break;
    end  
    if(sum((sum(U,2)>0))<=1)
        break;
    end

    Cond_old = Cond;
    %Cond
    
    % Record the number of clusters ----------------
    [~, idx] = max(U);
    uL = unique(idx);
    numclust(iter) = length(uL);
'';
end
end
clearvars -except E_ P_ U_ Cond_ numclust iter_breaks;
%% end of function



function [E] = initialize_E(X,U,M)
[~,L] = max(U);
E = cell(1,size(U,1));
for i=1:size(U,1)
    XX = X(L == i, :);
%    E{i} = XX(1:M,:);
    
    m = sum(L==i);
    if(m > M)
        [~, E{i}] = nnmf(XX,M);
%        r = unique(randi(m,1,50));
%        E{i} = XX(r(1:M),:);
    elseif(m == M)
        E{i} = XX;
    elseif(m < M)
        %E{i} = zeros(M,size(X,2));
        [~, E{i}] = nnmf(XX,m);
    end
end

