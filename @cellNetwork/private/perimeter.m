%% perimeter
% Compute of the perimeter of a given cell.
%%

%% Syntax   
% a = perimeter(V,A);
%
%% Description 
% Calculate the perimeter of a cell.
%
%% Inputs
% * V - a 2-by-2 matrix with the starting and ending points.
% * A - a variable containing the curvature angle of the edge with the
% cordlength.
%
%% Outputs
% * P - the perimeter of the cell.
%
%% Examples
%  >> a = perimeter([0 0;1 1;1 0],[0 0 0])
%
%% TODO
%
%% See also
% * 
%
%% Author  
% Sebastien Besson
% email address : jdumais@oeb.harvard.edu
% October 2008; Last revision: December 17, 2008

function perimeter = perimeter(V,A)

dV = circshift(V,[-1 0])-V;
L = sqrt(dV(:,1).^2+dV(:,2).^2);
perimeter = sum(abs(L.*A./sin(A)));

end