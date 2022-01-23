function [overlap,dissimilarity] = IDOA(X)


%Perform IDOA on input data 'X'
%Yogev Yonatan, yogev.yn@gmail.com
% 'X'- input data matrix, rows correspond to species,
% columns correspond to observations.
%
% 'overlap and 'dissimilarity' are the output overlap and
% dissimilarity matrices beween the observations, respectievly.
% The values in column i in 'overlap' matrix indicate the
% overlap between the i'th observation and all the other observations.
% 'dissimilarity' matrix values represent the correspondant dissimilarity
% (rJSD) between observations.

O = @(x,y) 0.5*(sum(x)+sum(y)); %Overlap function
X = X./sum(X);% normalization
m = size(X,2);% Number of observations

overlap = zeros(m-1,m);
dissimilarity = zeros(m-1,m);

KLD = @(x,y) sum(x(x>0).*log(x(x>0)./y(x>0)));
rJSD=@(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2));

for i = 1:m
    for j = 1+i:m
        Shared_species = find(X(:,i)>0 & X(:,j)>0);
        overlap(j-1,i) = O(X(Shared_species,i),X(Shared_species,j));
        overlap(i,j)=overlap(j-1,i);
        xhat = X(Shared_species,i)./sum(X(Shared_species,i));
        yhat =X(Shared_species,j)./sum(X(Shared_species,j));
        dissimilarity(j-1,i) = rJSD(xhat,yhat);
        dissimilarity(i,j)=dissimilarity(j-1,i);
    end
end






