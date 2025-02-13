function n = getNumberSpanningTrees(A)
%GETNUMBERSPANNINGTREES Get the number of spanning trees for adjacency matrix A
%   Returns the number of spanning trees for the connected undirected simple
%   graph described by the adjacency matrix A. The number of spanning trees
%   is determined using Kirchhoff's matrix tree theorem (cf., e.g., Russell
%   Merris, 'Laplacian Matrices of Graphs: A Survey', Linear Algebra and
%   its Applications, vol. 197-198, p. 143-176, 1994).

%   Authors:
%     Matthias Hotz <matthias.hotz@tum.de>
%
% Copyright (c) 2015, Matthias Hotz and Technische Universitaet Muenchen
% Covered by the 3-clause BSD License (see LICENSE file for details).
%--------------------------------------------------------------------------

% NOTE: We do not check for connectivity of the graph!

A = lib.ensureUndirectedSimpleGraph(A) ~= 0; % Check A and set weights to 1
A = A + A';

D = diag(sum(A, 2));    % Degree matrix
L = D - A;              % Laplacian matrix

% Matrix tree theorem:
%  The number of spanning trees is a cofactor of the Laplacian matrix
% (Note that rounding "removes" numerical inaccuracies...)
n = round(det(L(2:end, 2:end)));
    
end
