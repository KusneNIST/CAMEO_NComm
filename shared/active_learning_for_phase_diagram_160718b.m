function [query, certainty, risks, full_cluster_labels, fu] = active_learning_for_phase_diagram_160718b(S, XY, R, measured_indices, Ua, cluster_labels)

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

Nl = size(S,1)-size(XY,1);
N = size(S,1);
full_index = 1:size(S,1);
% active learning prediction
measured_indices = [1:Nl Nl+measured_indices]; 
measured_unmeasured_indices = [measured_indices setdiff(full_index, measured_indices)];
W = S(measured_unmeasured_indices, measured_unmeasured_indices);
fl = Ua;
U_blind = zeros(N,size(Ua,2));
U_blind(measured_indices,:) = Ua;
[fu, ~] = harmonic_function(W, fl); % 
[query, ~, risks] = active_learning(1, U_blind, full(S), measured_indices); % ~ = acc_ML
[guess_probability, best_guess] = max(fu,[],2);

certainty = ones(N,1);
full_cluster_labels = ones(N,1);
certainty(setdiff(full_index,measured_indices)) = guess_probability;
full_cluster_labels(measured_unmeasured_indices) = [cluster_labels' best_guess'];
full_cluster_labels = split_cluster_graph(S,full_cluster_labels,R);
