function [D d beta] = transform(type, D, d, V, U, beta, negIdx)

    [c n] = size(d);
    
    switch type
        case 'NE' 
            %get the index to the cluster and the point that caused the
            %negative distance
            [clusters, points]=ind2sub([c n],negIdx);
            uniqueClusters = unique(clusters)';
            tmp = zeros(c,n);
            
            for i = uniqueClusters
                k = points(clusters == i);
                tmp(i,k) = V(i,:)*V(i,:)' - 2*V(i,k) + 1;
            end
            
            deltaBeta = max(max((-2.*d(negIdx))./tmp(negIdx)));
            d(negIdx) = d(negIdx) + (deltaBeta/2).*tmp(negIdx);
            D = D + deltaBeta * (ones(n)-eye(n));
            beta = beta + deltaBeta;
            
        case 'BS'
            
        case 'SU'
            MST = graphminspantree(sparse(D));
            MST = MST + MST';
            
            for idx=1:length(j)

                [clust, z]=ind2sub([c n],j(idx));
                
                p = find_prototype(U(clust,:),D);
          
                if z ~= p
                    D(p,z) = su_dist(z,p,D,MST); 
                    D(z,p) = D(p,z);
                end

               [clust,z]=ind2sub([c n],j(idx));
               [~,p] = sort(V(clust,:),'descend');
               
               for pi=1:1 %length(p)
                if z ~= p(pi)
                    %fprintf('%d %d - Before %f\n',z,pi,D(z,pi));
                    D(z,p(pi)) = su_dist(z,p(pi),D,MST);
                    %fprintf('%d %d - After %f\n',z,pi,D(z,pi));
                    D(p(pi),z) = D(z,p(pi));
                end
               end

            end
            
            d(d<0) = min_d;
            
        case 'PF'
            
        case 'EF'
            
        case 'LF'
        
        
    end
end

