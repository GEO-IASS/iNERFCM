%clear MATLAB workspace
clear
close all

%load Iris dataset and compute the sup norm squared dissimilarity
a=normrnd(repmat([1,1],100,1), .1);
b=normrnd(repmat([1,5],100,1), 0.1);
c=normrnd(repmat([4,3],100,1), 0.1);
X = [a;b;c];

%{
a = randi([-100 100],3,2);
b = randi([-100 100],3,2);

X = round([a;b]);
X = [   -55   -47;
   -40    69;
   -79   -63;
    96   -40;
    -6   -37;
   -93    76];
%}
D = squareform(pdist(X,'chebychev')).^2;

%[ST,pred] = graphminspantree(sparse(D));
%view(biograph(ST,[],'ShowArrows','off','ShowWeights','on'));

n = size(D,1);

%compute the normalized dissimilarity image from D
f = figure('Visible','off');imagesc(D);colormap('gray');colorbar;
print(f, '-djpeg', 'Results/3Clouds/Images/Iris.jpg');

%D = D([1:5 60:65 120:125],[1:5 60:65 120:125])
transforms = {'SU','BS','PF','EP','LF','SU'};
        
%% iRFCM configurations/options (those are the default values)
options.fuzzifier        = 2;
options.epsilon          = 0.0001;
options.maxIter          = 100;
options.initType         = 2;

%set the number of clusters to 3
c= 3;

%% Since RFCM failed we need to run iRFCM
% loop for every delta
for i=1:1 %length(transforms)
    options.transform = transforms{i};
    out = inerfcm(D,c,options);
    
    %save the partition matrix for this delta
    U = out.U;
    dlmwrite(sprintf('Results/3Clouds/Partitions/U_%s(%d).csv',transforms{i},c),U, 'delimiter',',');

    %save the induced dissimilarity image for this delta
    %Ref. J. Huband and J. Bezdek, “VCV2– Visual cluster validity,” Comput. Intell. Res. Front., 2008.
    uu = 1 - ((U'*U)./max(max(U'*U)));
    f = figure('Visible','off');imagesc(uu);colormap('gray');caxis([0 1]);
    print(f, '-djpeg', sprintf('Results/3Clouds/Images/UU_%s(%d).jpg',transforms{i},c));
end