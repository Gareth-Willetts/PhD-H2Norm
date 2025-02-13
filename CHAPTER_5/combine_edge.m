function [src, dst, weight] = combine_edge(src, dst, weight)
    tbl = table(src', dst');
    [~, ia, ic] = unique(tbl, 'stable');
    
    for i = 1:length(weight)
        % If i does not agree with the ith entry of the uniqueness vector
        if i ~= ic(i)
            % Sum ith weight and the ith unique weight
            weight(ic(i)) = weight(i) + weight(ic(i));
            % Zero out the extra weight
            weight(i) = 0;
            break
        end
    end

    src = src(ia);
    dst = dst(ia);
    weight = weight(weight~=0);
end