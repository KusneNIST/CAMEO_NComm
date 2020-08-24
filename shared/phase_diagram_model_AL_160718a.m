function [U, E, P, cluster_labels] = phase_diagram_model_AL_160718a(X, R, S, T, experiment)
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

N = size(X,1);
num_clust = min([N experiment.num_clusters]);
    
W = SimGraph_Full(X(:,:)', experiment.sim_sigma, experiment.distance_type);
idx = SpectralClustering(W, num_clust, 2);
    
% form the cluster label matrix U
k = length(unique(idx));
U = zeros(size(idx,1),k);
for j=1:length(idx)
    U(j,idx(j)) = 1;
end

% Use graph cut to clean up the data
Ui = graph_cut_160402a(X,U',S,experiment.graphcut_weight);
k = size(Ui,1);
[~,idx2] = max(Ui);
experiment.spectral_clustering_graph_cut_cluster_number = k;
%

% package data
data.spectra = X - min(X(:));
data.composition = R; % composition
data.coordinates = []; % xy coordinates from conversion of composition
data.connectivity = S; % graph
data.spectra_cluster_labels = Ui; % cluster membership matrix
data.spectra_sampling_coordinates = T;  % can be 2theta, q-value, d-spacing, etc.

% Parameters for GraphCut    
parameters.graphcut.wgc         = experiment.graphcut_weight;  % 1E5  
parameters.graphcut.alpha       = experiment.graphcut_balance; % 0.75

% Parameters for PCOMMEND
parameters.PCOMMEND.alpha        = 0.0001;         % Regularization Parameter. Range of alpha should be between 0 and 1.
parameters.PCOMMEND.changeThresh = 1e-5;         % Stopping criteria for the objective function. 
parameters.PCOMMEND.M            = experiment.num_endmembers;            % Number of Endmembers per endmember set. 
parameters.PCOMMEND.iterationCap = experiment.max_iter;  %1500       % Maximum number of iterations.     
parameters.PCOMMEND.C            = experiment.num_clusters;            % Number of clusters or number of endmember sets. 
parameters.PCOMMEND.m            = 2;            % Fuzzifier for membership. 
parameters.PCOMMEND.EPS          = 0.0001;       % Small positive constant. 
parameters.PCOMMEND.GRENDEL_lower_bound_cutoff = experiment.GRENDEL_lower_bound_cutoff;
parameters.spectral_clustering_graph_cut_cluster_number = experiment.spectral_clustering_graph_cut_cluster_number;

[E,P,U_,Cond]=Phase_Mapping_160715a(data,parameters);

if(~isempty(Cond))
    [~, mCi] = min(Cond);
    U = U_{mCi}';
    E = E{mCi};
    P = P{mCi};
    [~, cluster_labels] = max(U, [], 2);
else
    cluster_labels = idx2';
    U = Ui';
end
clearvars -except U E P cluster_labels;