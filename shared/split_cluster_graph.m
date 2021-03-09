function [L] = split_cluster_graph(S, idx, C)
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

ui = unique(idx);
S = S > 0;
l = 1;
L = idx * 0;
% 
for i=1:length(ui)
    c = find(idx == ui(i));
    t = {};
    t{1} = c(1);
    cc{i} = t;
    for j=2:length(c)
        c_ = c(j);
        found = 0;
        m = cc{i};
        for k=1:length(m)
            if( sum(S(c_,cc{i}{k})) ) 
                cc{i}{k} = [cc{i}{k} c_];
                found = 1;
            end
        end
        if(~found); cc{i}{k+1} = c_; end;
    end
    
    t = cc{i};
    % merging groups
    join = 1;
    while(sum(join(:)))
        % identify groups to merge
        join = zeros(length(t),length(t));
        for j=1:length(t)
            for k=1:length(t)
                if k>j
                    s = S(t{j},t{k});
                    join(j,k) = sum(s(:)) > 0;
                end
            end
        end
        % merge 2 groups
        if( sum(join(:)))
            [m, n] = find(join);
            t{m(1)} = [t{m(1)} t{n(1)}];
            kp = 1:length(t);
            kp = exclude(kp,n(1));
            t = t(kp);
        end
    end
    cc{i} = t;
    
    m = cc{i};
    for j=1:length(m)
            L(cc{i}{j}) = l;
            l = l+1;
%        end
    end
end

ind = find(L == 0);

if(size(C,2) == 3)
    xy = tern2cart(C);
else
    xy = C;
end
D = squareform(pdist(xy));
