function A = ensureUndirectedSimpleGraph(A)
%ENSUREUNDIRECTEDSIMPLEGRAPH Ensure an undirected simple graph
%   This function reports an error if the adjacency matrix A does not
%   represent an undirected simple graph with 2 or more edges. If
%   successful, it returns a lower triangular adjacency matrix.
%
%   Note that this function requires Matlab R2014a or later due to the use
%   of 'issymmetric' and 'istril'.

%   Authors:
%     Matthias Hotz <matthias.hotz@tum.de>
%
% Copyright (c) 2015, Matthias Hotz and Technische Universitaet Muenchen
% Covered by the 3-clause BSD License (see LICENSE file for details).
%--------------------------------------------------------------------------

assert(ismatrix(A) && size(A, 1) == size(A, 2), ...
    'Adjacency matrix must be square.');
assert(all(diag(A) == 0), ...
    'Graph must not contain self-loops.');
assert(issymmetric(A) || istril(A), ...
    'Graph must be undirected (Adjacency matrix is symmetric or lower triangular).');
assert(~isscalar(A), ...
    'Graph must have two or more vertices.');

A = tril(A);
    
end
