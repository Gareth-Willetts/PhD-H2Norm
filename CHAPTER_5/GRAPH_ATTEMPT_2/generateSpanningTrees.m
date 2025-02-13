function [idx, src, dst] = generateSpanningTrees(A)
%GENERATESPANNINGTREES Generate all spanning trees for an undirected graph
%   This function generates all spanning trees for the undirected simple
%   graph that is described by the adjacency matrix A. It is an
%   implementation of 'Algorithm S' on p. 464 in The Art of Computer
%   Programming, Volume 4A (Combinatorial Algorithms, Part 1), by Donald E.
%   Knuth.
%
%   The output comprises an index matrix idx, which, for every spanning
%   tree, comprises the indices of the source and destination vertices. For
%   example, the following code outputs, for every spanning tree, a matrix
%   with the vertices of all its edges:
%
%    [idx, src, dst] = lib.generateSpanningTrees(A);
%    for n = 1:size(idx, 2)      % Iterate through all spanning trees
%        % Output source and destination vertex of all edges of spanning tree n:
%        [src(idx(:, n)), dst(idx(:, n))]
%    end

%   Authors:
%     Matthias Hotz <matthias.hotz@tum.de>
%
% Copyright (c) 2015, Matthias Hotz and Technische Universitaet Muenchen
% Covered by the 3-clause BSD License (see LICENSE file for details).
%--------------------------------------------------------------------------

% NOTE 1: Beforehand, I apologize for the rather tightly bounded beauty of
% the code below. I hope the tight bound gave rise to more efficiency.
% NOTE 2: This code does not include the improvements of 'Algorithm S' as
% described on p. 467 in TAOCP, vol. 4A, and the heurists for the bridge
% test, cf. p. 803 in TAOCP, vol. 4A, are omitted as well. For my purposes,
% they were not necessary, so I preferred to keep the code as slim as
% possible. However, you are welcome to add them ;)
% NOTE 3: When the number of spanning trees is large, returning the result
% as done here is very memory demanding. In this case you may, e.g., modify
% the code to use a callback function to process a newly found spanning
% tree in order to avoid storing them all. The section of the code, where
% the result is stored, is highlighted.

% Preparation: Generate data and define operations on it
[V, E, d, t, n, p, mate, src, dst, idx_map] = generateGraphData(A); % Local function: See bottom of file

function delete(a)
    n(p(a)) = n(a);
    p(n(a)) = p(a);
end

function undelete(a)
    p(n(a)) = a;
    n(p(a)) = a;
end

% S1: Initialization
skip_s2_to_s4 = false;  % Flag: Skips S2, S3, and S4 if true (used by S1)
i = 0;                  % Counter for spanning trees
total_num_st = lib.getNumberSpanningTrees(A);
idx = zeros(V - 1, total_num_st); % Edge index matrix for spanning trees

a = generateSpanningTree(V, t, n); % Local function: See bottom of file
s = zeros(V - 2, 1);    % Stack for temporarily removed arcs
l = zeros(V + 2*E, 1);  % Links to successors
level = 1;              % Current level (Called 'l' in TAOCP)
x = 0;                  % Arc node that was replaced in last ST to obtain current ST
                        % (x is not used in the result, but maybe we need it in the future ;) )

if V == 2
    v = 1;
    e = n(v);
    skip_s2_to_s4 = true;
end

% S2 to S9: Loop to repeat S2
while true
    
    if ~skip_s2_to_s4
        
        % S2: Enter level 'level'
        e = a(level + 1);   % Choose edge that should be part of the ST and
        u = t(e);           %  set u and v to the adjacent vertices
        v = t(mate(e));
        if d(u) > d(v)      % Ensure a higher (or equal) degree at the tail vertex
            [v, u] = deal(u, v);
            e = mate(e);
        end

        % S3: Contract edge 'e' (v -> u) and combine vertex u into vertex v
        k = d(u) + d(v);    % Degree of the combined vertex
        f = n(u);
        g = 0;
        while t(f) ~= 0     % Go through the arc node list of u
            if t(f) == v        % If the edge leads to v, remove it (not part of reduced graph)
                delete(f);
                delete(mate(f));
                k = k - 2;
                l(f) = g;       % Add link to recently removed arc g
                g = f;          % Set recently removed arc to f
            else
                t(mate(f)) = v; % Mate edge now leads to v (instead of u)
            end
            f = n(f);
        end
        l(e) = g;           % Remember stack of arcs removed for contraction of 'e'
        d(v) = k;
        n(p(u)) = n(v);     % Insert the reduced node arc list of u at
        p(n(v)) = p(u);     %  the beginning of the node arc list of v
        p(n(u)) = v;
        n(v) = n(u);
        a(level) = e;       % Edge e is part of the spanning tree (on this level)

        % S4: Move down one level
        level = level + 1;
        if level < V - 1    % Until the reduced graph has only 2 vertices:
            s(level) = 0;       % Clear 'removed arc'-stack of new level and...
            continue;           %  ...repeat S2
        end
        e = n(v);
        
    end % End: Skip S2 to S4

    % S5: Reduced graph comprises two vertices now, visit all its spanning trees
    while true
        a(V - 1) = e;           % Update edge on highest level of the spanning tree
        
        %----- NEW SPANNING TREE IN 'a': SAVE THE RESULT -----
        i = i + 1;              % Increase spanning tree counter
        idx(:, i) = idx_map(a); % Convert arc node indices to edge indices and save result
        %-----------------------------------------------------
        
        x = e;                  % Save the previously processed edge...
        e = n(e);               %  ...and take the next one
        if t(e) == 0
            break;              % We processed the entire arc list, so move on to S6
        end
    end

    % S6 to S9: Loop to repeat S6
    while true
        
        % S6: Move up one level
        level = level - 1;
        if level == 0
            assert(i == total_num_st, ['Something went wrong.\n' ...
                'Found %i spanning trees, but there should be %i.'], i, total_num_st);
            return;         % We are done!
        end
        e = a(level);
        u = t(e);
        v = t(mate(e));

        % S7: Revert contraction of edge 'e' (v -> u)
        n(v) = n(p(u));     % Remove the part of the arc node list of v
        p(n(v)) = v;        %  that was inherited from u during edge contraction
        n(p(u)) = u;        % ...and let both ends of the extracted part
        p(n(u)) = u;        %  point to the header node of u
        f = p(u);
        while t(f) ~= 0     % Go through the arc node list of u
            t(mate(f)) = u;     % Mate edge leads to u again (instead of v)
            f = p(f);
        end
        f = l(e);
        k = d(v);
        while f ~= 0        % Go through stack of arcs removed during contraction of 'e'
            k = k + 2;
            undelete(mate(f));  % ...and restore the removed edges
            undelete(f);
            f = l(f);
        end
        d(v) = k - d(u);
        
        % S8: If 'e' is not a bridge, remove it from the graph and repeat S2
        if ~isBridge(mate(e), u, v, V, t, n)
            x = e;
            l(e) = s(level);    % Put 'e' on top of the 'removed arc'-stack...
            s(level) = e;
            delete(e);          % ...and remove it from the graph
            delete(mate(e));
            d(u) = d(u) - 1;
            d(v) = d(v) - 1;
            break;              % Repeat S2 (by exiting the S6 loop) for reduced graph
        end
    
        % S9: Undo edge deletions of level 'level'
        e = s(level);
        while e ~= 0            % Go through the 'removed arc'-stack for this
            u = t(e);           %  level and recover the edges
            v = t(mate(e));
            d(u) = d(u) + 1;
            d(v) = d(v) + 1;
            undelete(mate(e));
            undelete(e);
            e = l(e);
        end
    end % S6
end % S2

end

function [V, E, d, t, n, p, mate, src, dst, idx_map] = generateGraphData(A)
%GENERATEGRAPHDATA Generate the graph data for Algorithm S
%   This function generates the data for the undirected simple graph
%   described by the adjacency matrix A as required by Algorithm S on
%   p. 464 in TAOCP, vol. 4A.
%
%   Note: The arc node index is 1-based here, but 0-based in TAOCP.

    A = lib.ensureUndirectedSimpleGraph(A) ~= 0; % Lower triangular & reset weights to 1
    [dst, src] = find(A);       % Directed edges: src(i)->dst(i) & dst(i)->src(i)
                                % (Note the sorting: src(i) < dst(i) & src(i) <= src(i+1) )
    V = size(A, 1);             % Number of vertices (in TAOCP: n)
    E = length(src);            % Number of edges

    d = full(sum(A, 2) + sum(A, 1)'); % Degree vector
    t = [src'; dst'];           % Tip vector
    t = [zeros(V, 1); t(:)];
    n = zeros(V + 2*E, 1);      % Next element in arc list
    p = zeros(V + 2*E, 1);      % Previous element in arc list
    for v = 1:V
        a = [V + 2*(find(dst == v) - 1) + 1; V + 2*find(src == v)];
        n([v; a]) = [a; v];
        p([a; v]) = [v; a];
    end
    mate = V + [2*(1:E); 2*(0:E-1) + 1]; % Maps from arc node index of an edge to
    mate = [zeros(V, 1); mate(:)];       % the arc node index of its "mate" edge
    idx_map = [1:E; 1:E];       % Maps from arc node index to edge index in src/dst
    idx_map = [zeros(V, 1); idx_map(:)];

    % FOR CODE ANALYSIS ONLY!
    %fprintf('\nd |'); fprintf(' %2d', d);
    %fprintf('\na |'); fprintf(' %2d', 1:(V+2*E));
    %fprintf('\nt |'); fprintf(' %2d', t);
    %fprintf('\nn |'); fprintf(' %2d', n);
    %fprintf('\np |'); fprintf(' %2d', p);
    %fprintf('\nm |'); fprintf(' %2d', mate);
    %fprintf('\n');
end

function a = generateSpanningTree(V, t, n)
%GENERATESPANNINGTREE Generates a spanning tree
%   This function implements the algorithm described in the solution of
%   exercise 94 on p. 802 in TAOCP, vol. 4A, to generate a spanning tree.
    
    a = zeros(V - 1, 1);    % Stores the arc pointers for the spanning tree
    b = zeros(V, 1);
    v = 1;                  % Start at (header node of) vertex 1
    w = 1;                  % Last visited vertex: Vertex 1
    b(1) = 1;               % Indicate vertex 1 was visited
    k = V - 1;              % We need to visit V - 1 more vertices

    while true
        e = n(v);           % Get first edge in vertex v's node arc list
        while t(e) ~= 0     % Continue until end of node arc list
            u = t(e);           % Get tip vertex of this edge
            if b(u) == 0        % If vertex u was not visited yet:
                b(u) = w;           % Indicate u was visited & save previously visited vertex
                w = u;              % Update the last visited vertex
                a(k) = e;           % Add this edge to the spanning tree
                k = k - 1;          % One vertex less on our 'to-visit-list'
                if k == 0           % All vertices were visited?
                    return;
                end
            end
            e = n(e);           % Process the next arc in v's node arc list...
        end
        if w == 1           % We were not able to visit all V vertices by following edges
            error('Graph must be connected.');
        end
        v = w;              % Move on to (header node of) vertex that was added last
        w = b(w);           % If w doesn't lead to unvisited vertices, continue with previous vertex
    end
end

function flag = isBridge(e_mate, u, v, V, t, n)
%ISBRIDGE Performs the bridge test for edge 'e' (v -> u), for which 'e_mate' is its mate
%   This function implements the algorithm described in the solution of
%   exercise 95 on p. 803 in TAOCP, vol. 4A, to check if 'e' is a bridge.

    b = zeros(V, 1);        % Stores the "vertex check list", which contains
                            % the vertices that are reachable from u but
                            % have not yet been checked if they lead to v
    w = u;                  % We are at vertex u
    b(w) = v;               % Indicate the end of the "vertex check list"

    while true
        f = n(u);
        while t(f) ~= 0     % Go through the arc node list of u
            v_prime = t(f);
            if b(v_prime) == 0  % If vertex v_prime was not visited yet:
                if v_prime ~= v     % Add v_prime to the "vertex check list"
                    b(v_prime) = v;
                    b(w) = v_prime;
                    w = v_prime;
                elseif f ~= e_mate  % Edge is not e_mate? We found a path from u to v!
                    flag = false;       % => Edge e/e_mate is not a bridge
                    return;
                end
            end
            f = n(f);
        end
        u = b(u);           % Process the next vertex of the "vertex check list"
        if u == v           % End of "vertex check list"?
            flag = true;        % ...then v is only reachable from u via via e_mate!
            return;
        end
    end
end
