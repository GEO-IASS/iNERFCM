%clear MATLAB workspace
clear
close all

%load animal mutation dataset
load Data/animal_mutation.csv;
D = animal_mutation;
n = size(D,1);

%compute the normalized dissimilarity image from D
D01 = D./max(D(:));
f = figure('Visible','off');imagesc(D.^2);colormap('gray');
print(f, '-djpeg', 'Results/Mutation/Images/animal_mutation.jpg');

%labels assigned to every delta
transforms = {'NE','BS','PF','EP','LF','SU'};
        
%% iRFCM configurations/options (those are the default values)
options.fuzzifier        = 2;
options.epsilon          = 0.0001;
options.maxIter          = 100;
options.initType         = 2;

%set the number of clusters to 4
c = 4;
    
%% loop for every delta
for i=1:1%length(transforms)
    options.transform = transforms{i};
    out = inerfcm(D.^2,c,options);

    %save the partition matrix for this delta
    U = out.U;
    dlmwrite(sprintf('Results/Mutation/Partitions/U-%s(%d).csv',transforms{i},c),U, 'delimiter',',');

    %save the induced dissimilarity image for this delta
    %Ref. J. Huband and J. Bezdek, “VCV2– Visual cluster validity,” Comput. Intell. Res. Front., 2008.
    uu = 1 - ((U'*U)./max(max(U'*U)));
    f = figure('Visible','off');imagesc(uu);colormap('gray');caxis([0 1]);
    print(f, '-djpeg', sprintf('Results/Mutation/Images/UU-%s(%d).jpg',transforms{i},c));
end