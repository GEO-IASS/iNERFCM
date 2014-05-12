function output = inerfcm(R, c, options)
%% 
%   Improved Relational Fuzzy c-Means (iRFCM) is an extension on top of the
%   Relational Fuzzy c-Means (RFCM) proposed in [1]. Since RFCM is the
%   relational dual of Fuzzy c-Means (FCM), it expects the input
%   relational matrix R to be Euclidean. Otherwise, RFCM can fail to
%   execute due encountering negative relational distances. iRFCM attempts
%   to solve this problem by Euclideanizing R first before clustering.
%
% Usage: output = irfcm(R,c,options)
%   options is a struct with the following default values:
%
%       fuzzifier        = 2;
%       epsilon          = 0.001;   
%       maxIter          = 100;     
%       initType         = 2;       
%       gamma            = 0;       
%       delta            = [];
%
%   Explanation of those fields is provided below
%
% output    - structure containing:
%               U: fuzzy partition
%               V: cluster centers/coefficients
%               terminationIter: the number of iterations at termination
%               maxIter: maximum number iterations allowed
%
%             If you decide to Euclideniate R before running iRFCM, then the
%             following information will be returned as well
%               kruskalStress: the stress value measured the tranformed R
%               eps: the value that minimized the different between R and D
%               lambda: the smallest constant to add to R to make it Euclidean
%               D: The Euclideaniated matrix
%
% R         - the relational (dissimilarity) data matrix of size n x n
% c         - number of clusters
% fuzzifier - fuzzifier, default 2
% epsilon   - convergence criteria, default 0.0001
% initType  - initialize relational cluster centers V
%               1 = random initialization
%               2 = randomly choose c rows from D
% maxIter   - the maximum number fo iterations, default 100
% delta     - delta is the matrix that is used to tranform R to Euclidean
%               (See ref. [2-3])
%             delta is of size n x n and can be for instance
%                  delta = 1-eye(n)  this is the default delta
%                  delta = subdominant ultrametric matrix
%                  delta = R
%                  delta = power of R
%                  delta = log2(1+D), see [3], page 485
%                  delta = parametric function
%              delta is expected to have an Euclidean representation,
%              otherwise an error will be thrown (see [2])
% gamma     - use it only if you know in advance that this the value that
%             will Euclideanize R
% 
% Refs:
%   [1] R. J. Hathaway and J. C. Bezdek, “Nerf c-means: Non-Euclidean 
%       relational fuzzy clustering,” Pattern Recognition, vol. 27, no. 3, pp. 429–437, Mar. 1994.
%   [2] J. Benasseni, M. B. Dosse, and S. Joly, “On a General Transformation Making a Dissimilarity 
%       Matrix Euclidean,” Journal of Classification, vol. 24, no. 1, pp. 33–51, Jun. 2007.
%   [3] J. Dattorro, Convex optimization and Euclidean distance geometry. 2005.

    %% iRFCM default values
    m = 2; epsilon = 0.0001;maxIter = 100;transformType='NE';
    
    %% Overwrite iRFCM options by the user defined options
    if nargin == 3 && isstruct(options)
        fields = fieldnames(options);
        for i=1:length(fields)
           switch fields{i}
               case 'fuzzifier', m = options.fuzzifier;
               case 'epsilon', epsilon = options.epsilon; 
               case 'initType', initType = options.initType; 
               case 'maxIter', maxIter = options.maxIter;
               case 'transform', transformType = options.transform;
           end
        end
    end
    
    %% Initialize variables
    D = R;n=size(D,1);d = zeros(c,n);bcount = 0;fcount = zeros(1,c);
    numIter=0;stepSize=epsilon;beta=0.0001; %U=Inf(c,n);
    
    %initialize relational cluster centers
    U = init_centers(initType, n, c, D);
    %U = [0.0570    0.2666    0.1265    0.2469         0    0.3030;
    %0.1932    0.1567    0.2594         0    0.0881    0.3026;
     %    0    0.2475    0.0106    0.4194    0.0442    0.2783];
    
    MST = graphminspantree(sparse(R));
    MST = MST + MST';
    changes = [];
    
    %% Begin the main loop:
    while  numIter < maxIter && stepSize >= epsilon
        U0 = U;
        
        V=U.^m;  
        V = V./(sum(V,2) * ones(1,n));
        
        U;
        V;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute the relational distances between "clusters" and points
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i=1:c
            d(i,:)=D*V(i,:)'-V(i,:)*D*V(i,:)'/2;
        end
        d;
        D;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Check for failure, are any of the d < 0?
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        negIdx = find(d(:) < 0)';
        if ~isempty(negIdx)
           fprintf('t=%d: found %d negative relational distances.\n',numIter, length(negIdx));
           
           %tranform the distance matrices here
           [D d beta, fcount, changes] = transform(transformType,R,D,d,V,U,beta,negIdx, fcount, MST, changes);
           bcount = bcount + 1;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Update the partition matrix U
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %First, compute U for only those k points where d > 0
        [~, k] = find(d > 0);
        
        d=d.^(1/(m-1));
        tmp = sum(1./d(:,k));
        U = zeros(c,n);
        U(:,k) = (1./d(:,k))./(ones(c,1)*tmp);
        
        %Second, for the points with d = 0
        %find the clusters and the points where d = 0
        [clusters, points] = find(d == 0);
        uniquePoints = unique(points)';
        
        for k = uniquePoints
            %some k might have a zero distance to more than one cluster
            idx = find(points == k);
            sub = sub2ind([c n],clusters(idx),points(idx));
            
            %The membership is 1/number of clusters to which k has d = 0
            U(sub) =  1/length(idx);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Update cluster prototypes V
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %V=U.^m;  
        %V = V./(sum(V,2) * ones(1,n));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Update the step size
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		stepSize=max(max(abs(U-U0)));
        
        numIter = numIter + 1;
    end
    
    %prepare output structure
    output = struct('U',U,...
                    'V',V,...
                    'terminationIter',numIter,...
                    'blockerCount',bcount,...
                    'D',D,...
                    'changes',changes);
                
    if nargin == 3,output.options = options;end
end