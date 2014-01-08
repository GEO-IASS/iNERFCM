function sdu = su_dist(i,j,D,MST)
%%
%
% Computes the Sub-dominant ultrametric matrix for a given dissimilarity
% matrix D based on its minimum spanning tree MST. If MST argument is not
% provided it will be computed assuming the MATLAB function graphminspantree
% exists. You need to have the Bioinformatics toolbox for that function.
%
% Usage sdu = subdominant_ultrametric(D,MST)
%
% sdu   - the computed Sub-dominant Ultrametric matrix
% D     - n x n dissimilarity matrix
% MST   - minimum spanning tree as returned by the MATLAB function 
%         graphminspantree

    %for every two nodes get the path between them
    [~, path, ~] = graphshortestpath(MST,i,j);
    pathLen = length(path) - 1;
    tmp = zeros(1,pathLen);
            
    %use the path to find the length of the edges
    for k=1:pathLen
        node = path(k:k+1);
        tmp(k) = D(node(1),node(2));
    end
            
    %find the longest edge, the largest distance
    sdu = max(tmp);
end