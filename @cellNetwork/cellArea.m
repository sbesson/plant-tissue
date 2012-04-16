%% cellArea
% Compute the oriented area of cells.
%%
%% Syntax   
% cellArea = cellArea(N,n)
%
%% Description
% This function calculates the oriented area of the cells of a cellNetwork
% object. The areas are oriented according to the right hand rule.
%
%% Inputs
% * N - a cellNetwork object
% * n - an array of cell indices
%
%% Outputs
% * cellArea - an array of the size of n containing the oriented areas of
% the corresponding cells.
%
%% Examples
% >> a =  cellArea(N,1);
% returns the area of the first cell of the N cellNetwork object. 
% >> a =  cellArea(N);
% returns the area of all the cells of the N cellNetwork object. 
%
%% See also 
% * cellCentroid
%
%% Author 
% Sebastien Besson.
% email address : sbesson@oeb.harvard.edu
% October 2007; Last revision: January 25, 2007

function [cellArea] = cellArea(N,n)

% Check the number of inputs
error(nargchk(1, 2, nargin));

% If no optional argument, loops over the cells
if (nargin == 1), n = 1:length(N.c); end;

% Initialize the output array
cellArea = zeros(length(n),1);

for c = 1:length(n)
    nc = n(c);
    
    %Extracts the cell data
    C = N.c{nc};
    
    % Extracts the vertices and angles
    V = zeros(length(C),2);
    A = zeros(length(C),1);
    
    index1 = logical(C>0);
    V(index1,:) = N.v(N.e(C(index1),1),:);
    A(index1) = N.e(C(index1),3); 
    
    index2 = logical(C<0);
    V(index2,:) = N.v(N.e(-C(index2),2),:);
    A(index2) = - N.e(-C(index2),3); 

    cellArea(c) = area(V,A);
end
end