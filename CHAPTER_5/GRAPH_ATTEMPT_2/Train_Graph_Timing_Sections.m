%% TO DO: generalise code!
clear all; close all; clc;
%% Create Graph
syms s mw ms mb kw cw x1 x2 kb ks
time_createG = [];
time_combineNodes = [];
time_admProd = [];
time_spanning = [];
admProdTba = [];
spanningTbaTime = [];

for j = 1:10000
    tic
    src = [1 1 2 2 3 4 1];
    dst = [3 4 3 3 4 5 5];
    Q1 = (ks/s) + x1;
    Q2 = (kb/s) + x2;
    weights_symbolic = [mw*s, mb*s, cw, kw/s, Q2, Q1, ms*s];
    
    %% Combine repeated edges - will replace this with generalised code later
    src = [1 1 2 3 4 1];
    dst = [3 4 3 4 5 5];
    
    weights_symbolic = [mw*s, mb*s, cw + kw/s, Q2, Q1, ms*s];
    time_createG(j) = toc;
    %% Combine vertices a and e, find admittance product
    tic
    [srcTae, dstTae, weights_symbolic_Tae] = combine_vertices(src, dst, weights_symbolic, 5,1);
    
    %% Combine repeated edges
    [srcTae, dstTae, weights_symbolic_Tae] = combine_edge(srcTae, dstTae, weights_symbolic_Tae);
    time_combineNodes(j) = toc;

    %% Finding spanning trees
    tic
    adj = createAdj(srcTae, dstTae);
    numSpanning_Tae = getNumberSpanningTrees(adj);
    
    % Find all spanning trees
    [idx, src_spn, dst_spn] = generateSpanningTrees(adj);
    % Create empty vector to store the source and destination pairs for each
    % spanning tree
    span_Tae = [];
    
    for n = 1:size(idx, 2)      % Iterate through all spanning trees
    % Store source and destination vertex of all edges of spanning tree n:
        span_Tae = [span_Tae;src_spn(idx(:, n)), dst_spn(idx(:, n))];
    end
    time_spanning(j) = toc;
    %% Combine vertices b and a
    [srcTba, dstTba, weights_symbolic_Tba] = combine_vertices(src, dst, weights_symbolic, 1,2);
    
    %% Combine repeated edges
    [srcTba, dstTba, weights_symbolic_Tba] = combine_edge(srcTba, dstTba, weights_symbolic_Tba);
    srcTba = srcTba - 1;
    dstTba = dstTba - 1;
    
    %% Finding spanning trees
    tic
    adj = createAdj(srcTba, dstTba);
    numSpanning_Tba = getNumberSpanningTrees(adj);
    
    % Find all spanning trees
    [idx, src_spn, dst_spn] = generateSpanningTrees(adj);
    % Create empty vector to store the source and destination pairs for each
    % spanning tree
    span_Tba = [];
    
    for n = 1:size(idx, 2)      % Iterate through all spanning trees
    % Store source and destination vertex of all edges of spanning tree n:
        span_Tba = [span_Tba;src_spn(idx(:, n)), dst_spn(idx(:, n))];
    end
    spanningTbaTime(j) = toc;
    
    %% Combine vertices b and e
    [srcTbe, dstTbe, weights_symbolic_Tbe] = combine_vertices(src, dst, weights_symbolic, 5,2);
    
    %% Combine repeated edges
    [srcTbe, dstTbe, weights_symbolic_Tbe] = combine_edge(srcTbe, dstTbe, weights_symbolic_Tbe);
    
    %% Finding spanning trees
    adj = createAdj(srcTbe, dstTbe);
    numSpanning_Tbe = getNumberSpanningTrees(adj);
    
    % Find all spanning trees
    [idx, src_spn, dst_spn] = generateSpanningTrees(adj);
    % Create empty vector to store the source and destination pairs for each
    % spanning tree
    span_Tbe = [];
    
    for n = 1:size(idx, 2)      % Iterate through all spanning trees
    % Store source and destination vertex of all edges of spanning tree n:
        span_Tbe = [span_Tbe;src_spn(idx(:, n)), dst_spn(idx(:, n))];
    end
    
    %% Undo transformation
    % span_Tba = span_Tba + 1; %undoing previous transformation
    % weights_symbolic_Tba = [0 weights_symbolic_Tba];
    
    %% Outline all edges and their corresponding weights
    % allCombinations = [1 2; 1 3; 1 4; 1 5; 2 3; 2 4; 2 5; 3 4; 3 5; 4 5];
    % weights = [0; mw*s; mb*s; ms*s; cw+kw/s; 0; 0; Q2; 0; Q1];
    edges_Tae = [srcTae' dstTae'];
    edges_Tba = [srcTba' dstTba'];
    edges_Tbe = [srcTbe' dstTbe'];
    
    %% Find admittance products in Tae and Tba
    tic
    sizeSpan_Tae = size(span_Tae,1) / numSpanning_Tae;
    for i = 1:sizeSpan_Tae:size(span_Tae,1)
        tree = span_Tae(i:i+sizeSpan_Tae-1,:); % find spanning tree
        [~, idx] = ismember(tree, edges_Tae, 'rows');
        if i == 1
            admittanceProds_Tae = prod(weights_symbolic_Tae(idx(idx~=0)));
        else
            admittanceProds_Tae = cat(3,admittanceProds_Tae, prod(weights_symbolic_Tae(idx(idx~=0))));
        end
    end
    
    tic
    sizeSpan_Tba = size(span_Tba,1) / numSpanning_Tba;
    for i = 1:numSpanning_Tba
        tree = span_Tba(i*sizeSpan_Tba - sizeSpan_Tba + 1:i*sizeSpan_Tba,:); % find spanning tree
        [~, idx] = ismember(tree, edges_Tba, 'rows');
        if i == 1
            admittanceProds_Tba = prod(weights_symbolic_Tba(idx(idx~=0)));
        else
            admittanceProds_Tba = cat(3,admittanceProds_Tba, prod(weights_symbolic_Tba(idx(idx~=0))));
        end
    end
    admProdTba(j) = toc;
    
    sizeSpan_Tbe = size(span_Tbe,1) / numSpanning_Tbe;
    for i = 1:numSpanning_Tbe
        tree = span_Tbe(i*sizeSpan_Tbe - sizeSpan_Tbe + 1:i*sizeSpan_Tbe,:); % find spanning tree
        [~, idx] = ismember(tree, edges_Tbe, 'rows');
        if i == 1
            admittanceProds_Tbe = prod(weights_symbolic_Tbe(idx(idx~=0)));
        else
            admittanceProds_Tbe = cat(3,admittanceProds_Tbe, prod(weights_symbolic_Tbe(idx(idx~=0))));
        end
    end
    
    
    %% Test
    tf = s*(sum(admittanceProds_Tae,3) + sum(admittanceProds_Tba,3) - sum(admittanceProds_Tbe,3))/(2*sum(admittanceProds_Tba,3));
    [n1,d1] = numden(tf);
    time_admProd(j) = toc;
end

%%
mean(time_createG(5:end))
mean(time_combineNodes(5:end))
mean(time_spanning(5:end))
mean(time_admProd(5:end))

%%
mean(admProdTba(5:end))
mean(spanningTbaTime(5:end))