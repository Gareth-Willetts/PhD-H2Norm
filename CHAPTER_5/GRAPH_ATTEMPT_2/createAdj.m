function adj = createAdj(src, dst)
    maxSrc = max(src); maxDst = max(dst);
    if maxSrc > maxDst
        adj = zeros(maxSrc, maxSrc);
    else
        adj = zeros(maxDst, maxDst);
    end

    for i = 1:length(src)
        adj(src(i), dst(i)) = 1;
        adj(dst(i), src(i)) = 1;
    end
end