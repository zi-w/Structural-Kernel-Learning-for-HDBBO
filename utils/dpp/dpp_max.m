function [C] = dpp_max(L, k)

C = [];

for i = 1:k
    [~, idx_max] = max(diag(L));
    
    if i < k
        beta = L(idx_max,:);
        L = L - beta' * beta / beta(idx_max);
    end

    C = [C,idx_max];
end

C = C(randperm(length(C)));

return