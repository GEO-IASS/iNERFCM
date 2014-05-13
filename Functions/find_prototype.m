function p = find_prototype(U,D)
    tmp = D*U';
    [~,p] = min(tmp);
end