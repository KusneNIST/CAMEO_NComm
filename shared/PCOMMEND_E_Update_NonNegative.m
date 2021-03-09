function E=PCOMMEND_E_Update_NonNegative(U,m,alpha,P,X,EPS)
% Modified by: A. Gilad Kusne, NIST, aaron.kusne@nist.gov
% Release 8/01/2020

%% This function updates the endmembers matrices (one per cluster)
% Input:
%   - U:      Fuzzy membership matrix (CxN). 
%   - m:      Fuzzifier.
%   - alpha:  Regularization Pprameter to trade off between the RSS and V
%             terms
%   - P:      Cell of C abundance matrices. One NxM matrix per cluster.
%   - X:      Pixel points (NxD matrix).
%
% Output:
%   - E:      Cell of C endmembers matrices. One MxD matrix per cluster
%
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
% THIS SOFTWARE IS PROVIDED BY THE UNIVERSITY OF FLORIDA AND
% CONTRIBUTORS ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES,
% INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED.  IN NO EVENT SHALL THE UNIVERSITY OR CONTRIBUTORS
% BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
% LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES,
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
% HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
% OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


% MODIFICATION COPYRIGHT:
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
%%

C    = length(P);
E    = cell(1,C);

for i=1:C
    E{i} = getEndmembers(U(i,:),m,P{i},X,alpha,EPS);
end
end
%% end of function

function E = getEndmembers(U,m,P,X,alpha,EPS)

D = size(X,2);
M = size(P,2);
N    = size(X,1);
DP   = EPS*eye(M,M);
Lmda = N*alpha/((M-1)*(1-alpha));
Z    = 2*Lmda*(eye(M,M)-(1/M)*ones(M,M));

try_E = 2;
if(M > 1)
    % Compute E.
    Y = (repmat((U.^m)', [1 M]).*P)';
    while try_E
    try
        E = pinv((Y*P+DP) + Z)*Y*X;   % Equation 6 in PCOMMEND paper
        try_E = 0;
    catch
        try_E = try_E - 1;
    end
    end
    % If some of the variables are negative, recurse.
    Ng = E' < 0;
    if(sum(Ng(:)) > 0)
        Ngu = unique(Ng,'rows','first');
        for i=1:size(Ngu,1)
            if(sum(Ngu(i,:)) > 0)
                vLocs = find(1 - Ngu(i,:)); % identify the positive variables in the set.
                % identify which variable (rows in Ng) have the same
                % pattern as Ngu(i,:).
                Ngui = repmat(Ngu(i,:),D,1);
                inds = all(Ng == Ngui,2); % boolean for rows with same patterns.
                % get new E for just the variables of interest.
                [Etemp] = getEndmembers(U,m,P(:,vLocs),X(:,inds),alpha,EPS); 
                Etemp2 = zeros(M,sum(inds));

                if(prod(size(vLocs))==0)
                end
                Etemp2(vLocs,:) = Etemp;
                E(:,inds) = Etemp2;
            end
        end
    end
else
    E = mean(X);  % if only one endmember for the set X, then set equal to first NNMF vector.
    %[~,E] = nnmf(X,1); % H is 1xD nonnegative factor.
end
end
