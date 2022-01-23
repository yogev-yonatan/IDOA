function [overlap,dissimilarity] = IDOA(x,X)

% Perform IDOA between observation 'x' and an observations' cohort 'X'.
% Yogev Yonatan, yogev.yn@gmail.com
% 'X'- Reference cohort, rows correspond to species,
% columns correspond to observations.
%
% 'overlap and 'dissimilarity' are the output overlap and
% dissimilarity values beween 'x' and the observations in 'X' , respectievly.
% The i'th value 'overlap' indicates the
% overlap between 'x' and the i'th observation in 'X'.
% The value in 'dissimilarity' represent the correspondant dissimilarity
% (rJSD) between 'x' and i'th observation in 'X'.


O = @(x,y) 0.5*(sum(x)+sum(y));% Overlap function
% normalizations:
X = X./sum(X);
x=x./sum(x);
%
m = size(X,2); %Number of observations in 'X'
%
overlap = nan(m,1);
dissimilarity = nan(m,1);


KLD = @(x,y) sum(x(x>0).*log(x(x>0)./y(x>0)));
rJSD=@(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2));

for i = 1:m    
    Shared_species = find(x>0 & X(:,i)>0);
    overlap(i) = O(X(Shared_species,i),x(Shared_species));
    xhat = X(Shared_species,i)./sum(X(Shared_species,i));
    yhat =x(Shared_species)./sum(x(Shared_species));
    dissimilarity(i) = rJSD(xhat,yhat);
end



