function [clust_perf_mp] = clustering_performance_measures_170930a(tl, el)

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

% assumed that tl is a row vector and el is a row matrix with each row corresponding to a
% clustering instance.
n = size(tl,2);
N = size(el,1);
tl_ = tl;
ul = unique(tl_);
for i=1:length(ul)
   tl(tl_ == ul(i)) = i; 
end

el_ = el;
ul = unique(el_);
for i=1:length(ul)
   el(el_ == ul(i)) = i; 
end

for k = 1:N
    % clustering performance
    new_labels_hat = max_precision_based_evaluation(tl, el(k,:));
    clust_perf_mp(k) = compute_perf(tl, new_labels_hat);%, stats.new_confusion_matrix);
end

function [perf] = compute_perf(tl, el)
n = size(tl,2);
td = squareform(pdist(tl'));
tsim = td == 0;
tdis = td ~= 0;

ed = squareform(pdist(el'));
esim = ed == 0;
edis = ed ~= 0;
TP_ = tsim & esim & ~eye(n,n);
FP_ = tsim & edis & ~eye(n,n);
FN_ = tdis & esim & ~eye(n,n);
TN_ = tdis & edis & ~eye(n,n);

TP = sum(TP_(:)) / 2;
FP = sum(FP_(:)) / 2;
FN = sum(FN_(:)) / 2;
TN = sum(TN_(:)) / 2;

% Fowlkes-Mallows Index
perf.fowlkes_mallows_index = TP / sqrt( (TP + FP) * (TP + FN) );

function [new_labels_hat] = max_precision_based_evaluation(labels, labels_hat)
[confusion_matrix,~] = confusionmat(labels,labels_hat);
[~, max_actual_classes] = max(confusion_matrix);
TPandFP = sum(confusion_matrix);

new_confusion_matrix = zeros( size(confusion_matrix,1) , size(confusion_matrix,1) );
list_of_combined_labels = cell(1,length(unique(labels)));
uL = unique([labels labels_hat]);
for i=1:length( max_actual_classes )
    best_matching_label = max_actual_classes(i);
    if(TPandFP(i) ~= 0)
        list_of_combined_labels{best_matching_label} = [list_of_combined_labels{best_matching_label} i];
    end
    new_confusion_matrix(:,best_matching_label) = new_confusion_matrix(:,best_matching_label) + confusion_matrix(:,i);
end

labels2drop = sum(new_confusion_matrix,1) == 0 & sum(new_confusion_matrix,2)' == 0;
new_confusion_matrix = new_confusion_matrix(~labels2drop,~labels2drop);

% new labels
list_of_combined_labels = list_of_combined_labels(:,sum(new_confusion_matrix) ~= 0);
new_labels_hat = zeros(size(labels_hat));
for i=1:size(list_of_combined_labels,1)
    for j=1:size(list_of_combined_labels,2)
        labels_to_fix = list_of_combined_labels{i,j};
        idx = false(size(new_labels_hat));
        for k=1:length(labels_to_fix)
            idx = idx | (labels_hat == labels_to_fix(k));
        end
        new_labels_hat(idx) = j;
    end
end
