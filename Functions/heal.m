function [D, d, changes] = transform(R, D, d, U, negIdx, MST, changes)

    [c, n] = size(U); 
    [clusters, points]=ind2sub([c n],negIdx);
    uniqueClusters = unique(clusters);
            
    for clust = uniqueClusters        
        p = find_prototype(U(clust,:),D);
        cp = find(clusters == clust);    
        k = points(cp);
            
        for z = k
            h = str2num(sprintf('%d%d%d%d',p,z,z,p));
                    
            if z ~= p && isempty(find(changes == h, 1))
                D(p,z) = su_dist(z,p,R,MST); 
                D(z,p) = D(p,z);
                changes = [changes;h];
            end
        end
                
    end
    
    d(negIdx) = d(negIdx) + -1*min(d(negIdx));
end

