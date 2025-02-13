function [tempsrc, tempdst, weight_new] = combine_vertices(src, dst, weight, nodechange, newnode)
    % Change nodechange -> newnode for T_(a,e)
    tempsrc = src;
    tempsrc(tempsrc==nodechange)=newnode;
    tempdst = dst;
    tempdst(tempdst==nodechange)=newnode;
    
    % Iterate through s (and t, as they are same length)
    for i = 1:length(tempsrc)
        % If both source and destination nodes are the new node
        if tempsrc(i) == newnode && tempdst(i) == newnode
            % Remove nodes and corresponding weight
            tempsrc = tempsrc(setdiff(1:end,i));
            tempdst = tempdst(setdiff(1:end,i));
            weight_new = weight(setdiff(1:end,i));
            % Break, as should only have 1 occurrence of this
            break
        else
            weight_new = weight;
        end
    end
    
    % Swap indices between src and dst so that lower numbers are always in src
    % This means that duplicate edges are spotted easier
    for i = 1:length(tempsrc)
        if tempsrc(i) > tempdst(i)
            temp = tempsrc(i);
            tempsrc(i) = tempdst(i);
            tempdst(i) = temp;
        end
    end
end