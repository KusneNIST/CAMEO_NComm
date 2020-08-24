function W = SimGraph_Full(M, sigma, metric)
% SIMGRAPH_FULL Returns full similarity graph
%   Returns adjacency matrix for a full similarity graph where
%   a Gaussian similarity function with parameter sigma is
%   applied.
%
%   'M' - A d-by-n matrix containing n d-dimensional data points
%   'sigma' - Parameter for Gaussian similarity function
%
%   Author: Ingo Buerk
%   Year  : 2011/2012
%   Bachelor Thesis

% Compute distance matrix
if(size(M,2) > 1)
    W = squareform(pdist(M',metric));
else
    W = 0;
end
% Apply Gaussian similarity function
W = simGaussian(W, sigma);

end