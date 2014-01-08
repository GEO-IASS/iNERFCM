function [D d beta] = transform(type, D, d, V, beta, j)

    [c n] = size(d);
    
    switch type
        case 'NE' %NERFCM approach
            min_d=1.0e-10;
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
               for k=1:n
                   [~, z]=ind2sub([c n],j(idx));
                   
                   if z ~= k
                        D(z,k) = su_dist(z,k,D,MST); 
                        D(z,k) = D(k,z);
                   end
               end
            end
            
        case 'PF'
            
        case 'EF'
            
        case 'LF'
        
        
    end
end

