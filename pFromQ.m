function p = pFromQ(Q)
p = transpose(null(transpose(Q)));

if size(p,1) > 1
    p = zeros(1,length(Q));
elseif size(p,1) == 1
    p = p./sum(p,2);
end


end