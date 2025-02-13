%%%%%%%%%%%%%%%%%% TRAIN SUSPENSION GRAPH %%%%%%%%%%%%%%%%%%%%%%%%%%%
%close all; clear all; clc
%%%%%%%%%%%%%%%%%% REFERENCES             %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matthias Hotz (2023). generateSpanningTrees(A) (https://www.mathworks.com/matlabcentral/fileexchange/53787-generatespanningtrees-a), MATLAB Central File Exchange. Retrieved October 9, 2023.
%% Create graph
syms s mw ms mb kw cw x1 x2 kb ks
% source and destination nodes
time = [];
time_make_graph = [];
time_spanning = [];
for i = 1:1
    tic
    src = [1 1 2 2 3 4 1];
    dst = [3 4 3 3 4 5 5];
    % arbitrary weights, so the adjacency matrix can be created later
    weights = [1, 1, 1, 1, 1, 1, 1];
    Q1 = (ks/s) + x1;
    Q2 = (kb/s) + x2;
    weights_symbolic = [mw*s, mb*s, cw, kw/s, Q2, Q1, ms*s];
    names = {'a' 'b' 'c' 'd' 'e'};
    % Create graph and plot
    G = graph(src,dst,weights,names);
    
    %% Combine repeated edges
    [src, dst, weights_symbolic_original] = combine_edge(src, dst, weights_symbolic);
    
    %% Combine vertices a and e, find admittance product
    [tempsrc, tempdst, weights_symbolic] = combine_vertices(src, dst, weights_symbolic_original, 5,1);
    
    %% Combine repeated edges
    [tempsrc, tempdst, weights_symbolic] = combine_edge(tempsrc, tempdst, weights_symbolic);
    time_make_graph(i) = toc;
    
    %% Finding number of spanning trees
    tic
    weights = ones(1,length(weights_symbolic));
    names = {'a' 'b' 'c' 'd'};
    G = graph(tempsrc,tempdst,weights,names);
    % Create adjacency matrix
    A = adjacency(G);
    % Store number of spanning trees
    numSpanning = getNumberSpanningTrees(A);
    
    % Find all spanning trees
    [idx, src_spn, dst_spn] = generateSpanningTrees(A);
    % Create empty vector to store the source and destination pairs for each
    % spanning tree
    answer = [];
    
    for n = 1:size(idx, 2)      % Iterate through all spanning trees
    % Store source and destination vertex of all edges of spanning tree n:
        answer = [answer;src_spn(idx(:, n)), dst_spn(idx(:, n))];
    end
    time_spanning(i) = toc;
    tic
    %% Compute the sum of products of each spanning tree
    T_ae = sum_admittance_products(answer, numSpanning, tempsrc, tempdst, weights_symbolic);
    
    %% Repeat for 2 -> 1
    %% Combine vertices b and a, find admittance product
    [tempsrc, tempdst, weights_symbolic] = combine_vertices(src, dst, weights_symbolic_original, 1, 2);
    
    %% Combine repeated edges
    [tempsrc, tempdst, weights_symbolic] = combine_edge(tempsrc, tempdst, weights_symbolic);
    % Have to subtract 1, as MATLAB requires one of the nodes to be
    % labelled as 1
    tempsrc = tempsrc - 1;
    tempdst = tempdst - 1;
    
    %% Finding number of spanning trees
    weights = ones(1,length(weights_symbolic));
    names = {'b' 'c' 'd' 'e'};
    G = graph(tempsrc,tempdst,weights,names);
    % Create adjacency matrix
    A = adjacency(G);
    % Store number of spanning trees
    numSpanning = getNumberSpanningTrees(A);
    
    % Find all spanning trees
    [idx, src_spn, dst_spn] = generateSpanningTrees(A);
    % Create empty vector to store the source and destination pairs for each
    % spanning tree
    answer = [];
    
    for n = 1:size(idx, 2)      % Iterate through all spanning trees
    % Store source and destination vertex of all edges of spanning tree n:
        answer = [answer;src_spn(idx(:, n)), dst_spn(idx(:, n))];
    end
    
    %% Compute the sum of products of each spanning tree
    T_ba = sum_admittance_products(answer, numSpanning, tempsrc, tempdst, weights_symbolic);
    
    %% Now repeat for T_(b,e)
    %% Combine vertices b and e, find admittance product
    [tempsrc, tempdst, weights_symbolic] = combine_vertices(src, dst, weights_symbolic_original, 5, 2);
    
    %% Combine repeated edges
    [tempsrc, tempdst, weights_symbolic] = combine_edge(tempsrc, tempdst, weights_symbolic);
    
    %% Finding number of spanning trees
    weights = ones(1,length(weights_symbolic));
    names = {'a' 'b' 'c' 'd'};
    G = graph(tempsrc,tempdst,weights,names);
    % Create adjacency matrix
    A = adjacency(G);
    % Store number of spanning trees
    numSpanning = getNumberSpanningTrees(A);
    
    % Find all spanning trees
    [idx, src_spn, dst_spn] = generateSpanningTrees(A);
    % Create empty vector to store the source and destination pairs for each
    % spanning tree
    answer = [];
    
    for n = 1:size(idx, 2)      % Iterate through all spanning trees
    % Store source and destination vertex of all edges of spanning tree n:
        answer = [answer;src_spn(idx(:, n)), dst_spn(idx(:, n))];
    end
    
    %% Compute the sum of products of each spanning tree
    T_be = sum_admittance_products(answer, numSpanning, tempsrc, tempdst, weights_symbolic);
    
    %% Combine for final answer
    final_answer = (1/2 * (T_ae + T_ba - T_be) * s) / T_ba;
    time(i) = toc;
end
time_full = time+time_spanning+time_make_graph;
mean(time_full(10:end))
mean(time_spanning(10:end))
mean(time_make_graph(10:end))