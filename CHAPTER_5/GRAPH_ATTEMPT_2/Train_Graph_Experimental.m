%% TO DO: generalise code!
clear all; close all; clc;
%% Create Graph
syms s mw ms mb kw cw x1 x2 kb ks
for i = 1:10000
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
    
    %% Combine vertices a and e, find admittance product
    [srcTae, dstTae, weights_symbolic_Tae] = combine_vertices(src, dst, weights_symbolic, 5,1);
    
    %% Combine repeated edges
    [srcTae, dstTae, weights_symbolic_Tae] = combine_edge(srcTae, dstTae, weights_symbolic_Tae);
    
    %% Finding spanning trees
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
    
    %%
    % for i = 1:numSpanning_Tba
    %     i*sizeSpan_Tba - sizeSpan_Tba + 1:i*sizeSpan_Tba
    % end
    
    
    
    %% Test
    tf = s*(sum(admittanceProds_Tae,3) + sum(admittanceProds_Tba,3) - sum(admittanceProds_Tbe,3))/(2*sum(admittanceProds_Tba,3));
    [n1,d1] = numden(tf);
    time(i) = toc;
end

%%
for i = 1:10000
    tic
    syms ms mb mw ks kb kw cw z s cb cs x y t;
    syms x1  x2 real
    x = [x1;x2];
    % Parameters for K1 and K2
    K1 = x1;
    K2 = x2;
    % Suspension admittances
    Q1 = (ks/s) + K1;
    Q2 = (kb/s) + K2;
    % State space representation of the system
    A2 = [ms*s^2+Q1*s, -s*Q1, 0; -s*Q1, s^2*mb+s*(Q1+Q2), -s*Q2; 0, -s*Q2, mw*s^2+Q2*s+cw*s+kw];
    B2 = [0; 0; cw*s+kw];
    % Extract transfer function
    TFT = [s, 0, 0]*inv(A2)*B2;
    TFTF = TFT(1,1);
    [n, d] = numden(TFTF);
    time_inbuilt(i) = toc;
end

mean(time(5:end))
mean(time_inbuilt(5:end))
%% Checking against working code
% simplify(n-n1)
% simplify(d-d1)
%% Checking against working code
% simplify(T_ae - sum(admittanceProds_Tae,3))
% simplify(T_ba - sum(admittanceProds_Tba,3))
% simplify(T_be - sum(admittanceProds_Tbe,3))

%% Stage 1
% for i = 1:size(admittanceProds_Tae,3)
%     for j = 1:size(admittanceProds_Tba,3)
%         simplify(admittanceProds_Tae(:,:,i) - admittanceProds_Tba(:,:,j))
%     end
% end
% No admittance products in common - can't demonstrate it for this example!

%% Stage 2

%% Stage 3

%% Stage 4

%% Stage 5

% Move all of this up to where the spanning trees are calculated -> then
% can work out the stages!


% %% Stage 0 - compare spanning trees to all possible edges and find number of occurrences
% %TaeANDTba = [];
% %R = [1 2 3 4];
% %allCombinations = combvec(R,R)';
% % src = [1 1 2 3 4 1];
% % dst = [3 4 3 4 5 5];
% % 
% % weights_symbolic = [mw*s, mb*s, cw + kw/s, Q2, Q1, ms*s];
% span_Tba = span_Tba + 1; %undoing previous transformation
% allCombinations = [1 2; 1 3; 1 4; 1 5; 2 3; 2 4; 2 5; 3 4; 3 5; 4 5];
% weights = [0; mw*s; mb*s; ms*s; cw+kw/s; 0; 0; Q2; 0; Q1];
% span_Tae_result = [];
% span_Tba_result = [];
% span_Tbe_result = [];
% 
% for i = 1:length(allCombinations)
%     temp = allCombinations(i,:);
%     span_Tae_result = [span_Tae_result; nnz(all(span_Tae == temp,2))];
%     span_Tba_result = [span_Tba_result; nnz(all(span_Tba == temp,2))];
%     span_Tbe_result = [span_Tbe_result; nnz(all(span_Tbe == temp,2))];
% end
% 
% %% Stage 1 - find admittance products in Tae and Tba
% TaeANDTba = min(span_Tae_result, span_Tba_result);
% 
% %% Stage 2 - find admittance products in Taa and Tbe = 0!
% Stage2 = zeros(length(TaeANDTba),1);
% 
% %% Stage 3 - find admittance products in stage 1 AND stage 2
% TaeANDTbaANDTbe = min(TaeANDTba, Stage2);
% 
% %% Stage 4 - admittance products obtained by removing stage 3 from stage 1
% Stage4 = TaeANDTba - TaeANDTbaANDTbe;
% 
% %% Stage 5 - admittance products obtained by removing stage 3 from stage 2
% Stage5 = Stage2 - TaeANDTbaANDTbe;
% 
% %% Calculate TF
% result = Stage4 - Stage5;
% finalsum = 0;
% for i = 1:length(result)
%     finalsum = finalsum + result(i) * weights(i);
% end
% 
% finalsum = finalsum*s;
% 
% finaldiv = 0;
% for i = 1:length(span_Tba_result)
%     finaldiv = finaldiv + 1/(span_Tba_result(i)*weights(i));
% end
% 
% finalsum*finaldiv


%%
% for i = 1:numSpanning_Tae
%     tempTae = span_Tae(i:i+size(span_Tae,1)/numSpanning_Tae-1,:); %find a spanning tree
%     for j = 1:numSpanning_Tba
%         tempTba = span_Tba(j:j+size(span_Tba,1)/numSpanning_Tba-1,:);
% 
%     end
% end

%% Stage 3
% sz = [size(span_Tae,1)/numSpanning_Tae, 2];
% span_Tae_tbl = table('Size', sz, 'VariableTypes', ["double", "double"]);
% finalTae = [];
% finalTbe = [];
% for i = 1:numSpanning_Tae %iterate through Tae spanning trees
%     %display(span_Tae(i:i+size(span_Tae,1)/numSpanning_Tae-1,:))
%     tempTae = span_Tae(i:i+size(span_Tae,1)/numSpanning_Tae-1,:); %find a spanning tree
%     for j = 1:numSpanning_Tbe %iterate through Tbe spanning trees
%         tempTbe = span_Tbe(i:i+size(span_Tbe,1)/numSpanning_Tbe-1,:); %find a spanning tree
%         tempTae = tempTae(not(ismember(tempTae, tempTbe, 'rows')),:); %find common admittance products
%         tempTbe = tempTbe(not(ismember(tempTae, tempTbe, 'rows')),:);
%         finalTbe = [finalTbe; tempTbe]; % problem here - this is returning way too many admittance products
%         if isempty(tempTae) || j == numSpanning_Tbe %if all have been cancelled out -> skip to next iteration, don't need to check any more
%             finalTae = [finalTae; tempTae];
%             break
%         end
%     end
% end

% for i = 1:numSpanning_Tbe
%     %display(span_Tbe(i:i+size(span_Tbe,1)/numSpanning_Tbe-1,:))
%     tempTbe = span_Tbe(i:i+size(span_Tbe,1)/numSpanning_Tbe-1,:);
% end

