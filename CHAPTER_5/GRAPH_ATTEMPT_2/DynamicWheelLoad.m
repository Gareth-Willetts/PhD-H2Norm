%% TO DO: generalise code!
clear all; close all; clc;
% v_ab,bc = T_ac + T_bb - T_ab - T_bc / 2T_ab
% T_bb = 0, T_ab = T_ba which has already been calculated, so need to find
% v_ab,bc = T_ac - T_ab - T_bc / 2T_ab
% Need to calculate T_ac and T_bc
%% Create Graph
syms s mw ms mb kw cw x1 x2 kb ks
cb = 7.1E3;
% for i = 1:10000
%     tic
src = [1 1 2 2 3 4 1];
dst = [3 4 3 3 4 5 5];
% Simpler example
K1 = 21065*s/(s+3.6373);
K2 = cb;
% K1 and K2 1st
% K1 = 1.2888E5*s/(s+2.4745);
% K2 = (701.13*s + 1.7926E9)/(s + 3.3969E3);
% K1 and K2 2nd
% K1 = (95319*s^2 + 4.4954E7*s + 3.5941E10)/(s^2 + 522.7415*s + 2.2583E5);
% K2 = (5622.9*s^2 + 8.7680E8*s + 5.8055E12)/(s^2 + 4933.1*s + 2.2439E7);
% K1 and K2 3rd
% K1 = (113.45*s^3 + 437810*s^2 + 9981.5*s + 0.000023395)/(s^3 + 5051.2*s^2 + 525540*s + 165.5255); %x1;
% K2 = (3.2132E12*s^3 + 5.8967E18*s^2 + 5.5411E9*s + 5.4456E-3)/(s^3 + 5.9915E12*s^2 + 160.1536*s + 6.4593E-10);%x2;
% Suspension admittances
Q1 = (ks/s) + K1;
Q2 = (kb/s) + K2;
weights_symbolic = [mw*s, mb*s, cw, kw/s, Q2, Q1, ms*s];

%% Combine repeated edges - will replace this with generalised code later
src = [1 1 2 3 4 1];
dst = [3 4 3 4 5 5];

weights_symbolic = [mw*s, mb*s, cw + kw/s, Q2, Q1, ms*s];

%% Combine vertices a and c, find admittance product
[srcTac, dstTac, weights_symbolic_Tac] = combine_vertices(src, dst, weights_symbolic, 1,3);

%% Combine repeated edges
[srcTac, dstTac, weights_symbolic_Tac] = combine_edge(srcTac, dstTac, weights_symbolic_Tac);
srcTac = srcTac - 1;
dstTac = dstTac - 1;

%% Finding spanning trees
adj = createAdj(srcTac, dstTac);
numSpanning_Tac = getNumberSpanningTrees(adj);

% Find all spanning trees
[idx, src_spn, dst_spn] = generateSpanningTrees(adj);
% Create empty vector to store the source and destination pairs for each
% spanning tree
span_Tac = [];

for n = 1:size(idx, 2)      % Iterate through all spanning trees
% Store source and destination vertex of all edges of spanning tree n:
    span_Tac = [span_Tac;src_spn(idx(:, n)), dst_spn(idx(:, n))];
end

%% Combine vertices b and a
[srcTba, dstTba, weights_symbolic_Tba] = combine_vertices(src, dst, weights_symbolic, 1,2);

%% Combine repeated edges
[srcTba, dstTba, weights_symbolic_Tba] = combine_edge(srcTba, dstTba, weights_symbolic_Tba);
srcTba = srcTba - 1;
dstTba = dstTba - 1;

%% Finding spanning trees
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

%% Combine vertices b and c
[srcTbc, dstTbc, weights_symbolic_Tbc] = combine_vertices(src, dst, weights_symbolic, 2,3);

%% Combine repeated edges
[srcTbc, dstTbc, weights_symbolic_Tbc] = combine_edge(srcTbc, dstTbc, weights_symbolic_Tbc);
srcTbc(srcTbc > 2) = srcTbc(srcTbc > 2) - 1;
dstTbc(dstTbc > 2) = dstTbc(dstTbc > 2) - 1;

%% Finding spanning trees
adj = createAdj(srcTbc, dstTbc);
numSpanning_Tbc = getNumberSpanningTrees(adj);

% Find all spanning trees
[idx, src_spn, dst_spn] = generateSpanningTrees(adj);
% Create empty vector to store the source and destination pairs for each
% spanning tree
span_Tbc = [];

for n = 1:size(idx, 2)      % Iterate through all spanning trees
% Store source and destination vertex of all edges of spanning tree n:
    span_Tbc = [span_Tbc;src_spn(idx(:, n)), dst_spn(idx(:, n))];
end

%% Undo transformation
span_Tbc(span_Tbc > 2) = span_Tbc(span_Tbc > 2) + 1;
% weights_symbolic_Tba = [0 weights_symbolic_Tba];

%% Outline all edges and their corresponding weights
% allCombinations = [1 2; 1 3; 1 4; 1 5; 2 3; 2 4; 2 5; 3 4; 3 5; 4 5];
% weights = [0; mw*s; mb*s; ms*s; cw+kw/s; 0; 0; Q2; 0; Q1];
edges_Tae = [srcTac' dstTac'];
edges_Tba = [srcTba' dstTba'];
edges_Tbe = [srcTbc' dstTbc'];

%% Find admittance products in Tae and Tba
sizeSpan_Tac = size(span_Tac,1) / numSpanning_Tac;
for i = 1:sizeSpan_Tac:size(span_Tac,1)
    tree = span_Tac(i:i+sizeSpan_Tac-1,:); % find spanning tree
    [~, idx] = ismember(tree, edges_Tae, 'rows');
    if i == 1
        admittanceProds_Tac = prod(weights_symbolic_Tac(idx(idx~=0)));
    else
        admittanceProds_Tac = cat(3,admittanceProds_Tac, prod(weights_symbolic_Tac(idx(idx~=0))));
    end
end

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

sizeSpan_Tbc = size(span_Tbc,1) / numSpanning_Tbc;
for i = 1:numSpanning_Tbc
    tree = span_Tbc(i*sizeSpan_Tbc - sizeSpan_Tbc + 1:i*sizeSpan_Tbc,:); % find spanning tree
    [~, idx] = ismember(tree, edges_Tbe, 'rows');
    if i == 1
        admittanceProds_Tbc = prod(weights_symbolic_Tbc(idx(idx~=0)));
    else
        admittanceProds_Tbc = cat(3,admittanceProds_Tbc, prod(weights_symbolic_Tbc(idx(idx~=0))));
    end
end

%% Test
tfnc = (sum(admittanceProds_Tac,3) - sum(admittanceProds_Tba,3) - sum(admittanceProds_Tbc,3))/(2*sum(admittanceProds_Tba,3));
tfnc = tfnc*(200/(200+s))*(kw/s + cw);
[n,d] = numden(tfnc);

%run("D:\OneDrive - University of Exeter\PhD\THESIS_WORKING_VERSION\CODE\CHAPTER_6\Train_J3_3rd")
run("C:\Users\Gareth\OneDrive - University of Exeter\PhD\THESIS_WORKING_VERSION\CODE\CHAPTER_6\Train_J3_3rd.m")