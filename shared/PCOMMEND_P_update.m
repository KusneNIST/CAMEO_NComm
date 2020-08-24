function P=PCOMMEND_P_update(X,E,C,EPS)
%% This function updates the abundace matrices (one per cluster)
% Input:
%   - X:   Pixel points (NxD matrix).
%   - E:   Cell of C endmembers matrices. One (MxD) matrix per cluster.
%   - M:   Number of endmembers.
%   - EPS: Small positive constant.
% Output:
%   - P:  Cell of C abundance matrices. One (NxM) matrix per cluster
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
%%



for i=1:C
    P{i} = getProps(E{i}, X, EPS);
end
end

function P = getProps(E, X, EPS)
% if one of the proportions are negative, try re-determining proportions based
% on the E's that provided positive values.
% do this recursively until only one E is left and the p value must be 1.
% the sum-to-one constraint keeps from all 3 being negative.
M  = size(E,1);
DP = EPS*eye(M,M);
N  = size(X,1);

% http://engineers.missouri.edu/zarea/files/Piecewise-Convex-Multiple-Model-Endmember-Detection-and-Spectral-Unmixing.pdf
if(M > 1)
    % P update
    Y = inv((E*(E)')+ DP); % M x M, correlation matrix between each of M 'endmembers'
    Z = ones(1,M)*Y; % 1 x M
    % Lambda is the Lagrange multiplier to enforce sum-to-one constraint
    % and maybe also Karush-Kuhn-Tucker condition for nonnegativity
    % constraints on the proportion values.
    Lamda=(1-(Z*E*(X)'))/(Z*ones(M,1));  % Equation 8 in PCOMMEND paper
    % update equation for proportion values based on objective function in
    % paper.
    P=(Y*(E*(X)'+ (Lamda'*ones(M,1)')'))';   % Equation 7 in PCOMMEND paper
    Z = P < 0;
    
    % If P<0, re-calculate proportions
    if(sum(sum(Z)) > 0)
        ZZ = unique(Z, 'rows', 'first');  % get unique first rows from Z
        for i = 1:size(ZZ, 1)
            if(sum(ZZ(i,:)) > 0)
                eLocs = find(1 - ZZ(i,:));
                rZZi = repmat(ZZ(i,:),N,1);
                inds = all(Z == rZZi,2);  % Find nonzero array in Z and store it in inds
                [Ptemp] = getProps(E(eLocs,:), X(inds,:), EPS);
                Ptemp2 = zeros(sum(inds), size(ZZ,2));
                
                Ptemp2(:, eLocs) = Ptemp;
                P(inds,:) = Ptemp2;
            end
        end
    end
else
    P = ones(size(X,1), 1);
    %P = max(X*(E'),0);
end
end


%% end of function