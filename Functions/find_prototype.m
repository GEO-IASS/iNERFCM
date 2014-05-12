function [p, k] = find_prototype(U,D, fcount, c)
    tmp = D*U';
    [~,p] = min(tmp);
    %[~, idx] = sort(tmp);
    
    k = find(U>0.4);
    %length(k)
    %[~, idx] = sort(U);out
    
    %k = idx(1:100)';
    
    t = max( max(U)/(1+log10(fcount)) , min(1/c, max(U)) );
    k = find(U>=t);
    
    %{
    if fcount == 2
        
    length(k)
    k
    U
    max(U)
    (1+log(fcount))
    t
    end
    %}
end