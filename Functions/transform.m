function [D d beta] = transform(type, D, d, V, U, beta, j)

    [c n] = size(d);
    min_d=1.0e-10;
    
    switch type
        case 'NE' 
            d_adjustment=zeros(c,n);
            
            for i=1:c
				work = (V(i,:) * V(i,:)' +1) / 2;
				d_adjustment(i,:) = work - V(i,:); 
            end
            
			work = (min_d - d(j)) ./ d_adjustment(j);
			beta_adjustment = max(work);
			beta = beta + beta_adjustment;
			d = d + beta_adjustment * d_adjustment;
			D = D + beta_adjustment * (ones(n) - eye(n));
            
            % Third, adjust all d values to be at least as big as min_d:
            d(d<min_d) = min_d;
            
        case 'BS'
            
        case 'SU'
            MST = graphminspantree(sparse(D));
            MST = MST + MST';
            
            for idx=1:length(j)
<<<<<<< HEAD
                [clust, z]=ind2sub([c n],j(idx));
                
                p = find_prototype(U(clust,:),D);
          
                if z ~= p
                    D(p,z) = su_dist(z,p,D,MST); 
                    D(z,p) = D(p,z);
                end
=======
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
>>>>>>> e266796c6cd380dd8acb3ee3bc4ac9ddc1ae456e
            end
            
            d(d<min_d) = min_d;
            
        case 'PF'
            
        case 'EF'
            
        case 'LF'
        
        
    end
end

