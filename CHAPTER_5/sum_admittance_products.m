function T = sum_admittance_products(answer, numSpanning, tempsrc, tempdst, weights, symbolic)
    arguments
        answer
        numSpanning
        tempsrc
        tempdst
        weights
        symbolic (1,1) double = 1
    end
    syms s;
    
    %% Compute the sum of products of each spanning tree
    % Currently only for spanning tree with 4 nodes, due to the symbolic
    % weights - need to generalise this
    
    % Find the size of each spanning tree = size of answer vector / number of
    % trees - each tree will be the same size!
    sizeGraph = size(answer,1)/numSpanning;
    % Create vector for the source and destination vertices
    % used to make original graph
    graphConnections = [tempsrc;tempdst];
    if symbolic
        %T = 0;
        answer_sym = cell(size(answer,1),1);
        % Compare spanning trees to original vertices
        for n = 1:size(answer,1)
            %answer_sym = 1;
            for i = 1:size(graphConnections,2)
                % if they match, append corresponding edge weight to the spanning
                % tree matrix
                if isequal(answer(n,:)', graphConnections(:,i))
                    %answer_sym{mod(n, sizeGraph)+1} = weights(i);%poly2sym(weights(i),s);
                    answer_sym{n} = weights(i);
                    break
                end
            end
            %if mod(n,sizeGraph) == 1 && n ~= 1
            %    T = T + answer_sym; %prod(cellfun(@prod,answer_sym, 'UniformOutput',false));
            %    answer_sym = 0;
            %end
        end
        
        % To compute product of spanning tree, need to create start and end vectors
        % for each spanning tree
        %startSeq = 1:sizeGraph:(numSpanning)*sizeGraph;
        %endSeq = sizeGraph:sizeGraph:numSpanning*sizeGraph;
        T = 0;
        
        % Calculate product of each spanning tree and sum them
        % stuffToAdd = 1;
        % counter = 1;
        for i = 1:sizeGraph:size(answer,1)
            % stuffToAdd = stuffToAdd*answer_sym{i};
            % if counter ~= sizeGraph
            %     counter = counter+1;
            % else
            %     T = T + stuffToAdd;
            %     counter = 1;
            %     stuffToAdd = 1;
            % end
            T = T + answer_sym{i}*answer_sym{i+1}*answer_sym{i+2};
        end
    else
    % Compare spanning trees to original vertices
        for n = 1:size(answer,1)
            for i = 1:size(graphConnections,2)
                % if they match, append corresponding edge weight to the spanning
                % tree matrix
                if answer(n,1:2)' == graphConnections(:,i)
                    answer(n,3) = weights(i);
                    break
                end
            end
        end
        
        % To compute product of spanning tree, need to create start and end vectors
        % for each spanning tree
        startSeq = 1:sizeGraph:(numSpanning)*sizeGraph;
        endSeq = sizeGraph:sizeGraph:numSpanning*sizeGraph;
        T = 0;
        
        % Calculate product of each spanning tree and sum them
        for i = 1:length(startSeq)
            T = T + prod(answer(startSeq(i):endSeq(i),sizeGraph));
        end
    end
end