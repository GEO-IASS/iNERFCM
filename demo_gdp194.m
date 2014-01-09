%clear MATLAB workspace
clear
close all

%load Iris dataset and compute the sup norm squared dissimilarity
D = load('Data/GDP194_FMS.csv');
n = size(D,1);

%compute the normalized dissimilarity image from D
D01 = D./max(D(:));
f = figure('Visible','off');imagesc(D.^2);colormap('gray');colorbar;
print(f, '-djpeg', 'Results/GDP194/Images/GDP194.jpg');

transforms = {'SU','SU','BS','PF','EP','LF'};
                    
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
    out = inerfcm(D.^2,c,options);
    
    %save the partition matrix for this delta
    U = out.U;
    dlmwrite(sprintf('Results/GDP194/Partitions/U_%s(%d).csv',transforms{i},c),U, 'delimiter',',');

    %save the induced dissimilarity image for this delta
    %Ref. J. Huband and J. Bezdek, “VCV2– Visual cluster validity,” Comput. Intell. Res. Front., 2008.
    uu = 1 - ((U'*U)./max(max(U'*U)));
    f = figure('Visible','off');imagesc(uu);colormap('gray');caxis([0 1]);
    print(f, '-djpeg', sprintf('Results/GDP194/Images/UU_%s(%d).jpg',transforms{i},c));
end