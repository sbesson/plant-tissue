%% centroid
% Conpute of the center of mass of a given cell.
%%

%% Syntax   
% C = centroid(V,A);
%
%% Description 
% Calculate the coordinate of the centroid of a cell.
%
%% Inputs
% * V - a 2-by-2 matrix with the starting and ending points.
% * A - a variable containing the curvature angle of the edge with the
% cordlength.
%
%% Outputs
% * C - the coordinates of the cell centroid.
%
%% Examples
%  >> C = centroid([0 0;1 1;1 0],[0 0 0])
%
%% TODO
%
%% See also
% * 
%
%% Author  
% Sebastien Besson
% email address : sbesson@oeb.harvard.edu
% October 2007; Last revision: December 22, 2008

function C = angularCentroid(V,A)
T = circshift(V,[-1 0])-V;
L = sqrt(T(:,1).^2+T(:,2).^2);
N= circshift(T,[0 1]);
N(:,2) = -N(:,2);

aL = zeros(size(V,1),1);
h = zeros(size(V,1),1);

for i=1:size(V,1)
    findP(V,A,s)
end
% Compute area of circular segments
ind1 = A.^2<=0.00001;
aL(ind1) = L(ind1).*(1+1/6*A(ind1).^2);
h(ind1) = A(ind1)/6;

ind2 = A.^2>0.00001;
aL(ind2) = L(ind2).*A(ind2)./sin(A(ind2));
h(ind2) = 1./A(ind2) -1./tan(A(ind2));

P = V+T/2+ repmat(h,1,2).*N/2;
C = sum(repmat(aL,1,2).*P,1)/sum(aL);
end