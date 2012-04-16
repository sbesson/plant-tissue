%% area
% Compute of the area of a given cell.
%%

%% Syntax   
% a = area(V,A);
%
%% Description 
% Calculate the area of a cell.
%
%% Inputs
% * V - a 2-by-2 matrix with the starting and ending points.
% * A - a variable containing the curvature angle of the edge with the
% cordlength.
%
%% Outputs
% * a - the oriented area of the cell.
%
%% Examples
%  >> a = area([0 0;1 1;1 0],[0 0 0])
%
%% TODO
%
%% See also
% * 
%
%% Author  
% Jacques Dumais
% email address : jdumais@oeb.harvard.edu
% October 2007; Last revision: December 17, 2008

function area = area(V,A)

dV = circshift(V,[-1 0])-V;
L = sqrt(dV(:,1).^2+dV(:,2).^2);

ind1 = A.^2<=0.00001;
subArea(ind1) = L(ind1).^2.*A(ind1)/6;
ind2 = A.^2>0.00001;
subArea(ind2) = L(ind2).^2.*(A(ind2)-cos(A(ind2)).*sin(A(ind2)))./(2*sin(A(ind2))).^2;

polArea = 1/2*(V(:,1).*circshift(V(:,2),-1)-V(:,2).*circshift(V(:,1),-1));
area = sum(polArea) + sum(subArea);

end