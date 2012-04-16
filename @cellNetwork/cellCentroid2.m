%% cellCentroid
% Compute the centroid of cells.
%%

%% Syntax
% centroid = cellCentroid(N,n)
% 
%% Description
% Compute the centroid of a given array of cells 
%
%% Inputs
% * N - a cellNetwork object
% * n - an array of cell indices
%
%% Outputs
% * centroid - a 2-element vector containing the coordinates of the centroid of the cell
%
%% Example
% >> cC =  cellCentroid(N,1);
% returns the centroid of the first cell of the N cellNetwork object. 
% >> cC =  cellCentroid(N);
% returns the centroid of all the cells of the N cellNetwork object. 
%
%% See also
% * cellArea
%
%% Author
% Sébastien Besson
% email address : jdumais@oeb.harvard.edu
% November 2007; Last revision: January 25, 2007

function [cC P] = cellCentroid2(N,n)

% Check the number of inputs
error(nargchk(1, 2, nargin));

% If no optional argument, loops over the cells
if (nargin == 1), n = 1:length(N.c); end;

% Initialize the output array
cC = zeros(length(n),2);

for c = 1:length(n)
    %Extract the cell data
    C = N.c{n(c)};
    
    % Extract the vertices and angles
    V = zeros(length(C),2);
    A = zeros(length(C),1);
    
    index1 = logical(C>0);
    V(index1,:) = N.v(N.e(C(index1),1),:);
    A(index1) = N.e(C(index1),3);
    
    index2 = logical(C<0);
    V(index2,:) = N.v(N.e(-C(index2),2),:);
    A(index2) = - N.e(-C(index2),3);

    % Compute the centroid of the cell
    [cC(c,:) P(:,:,c)] = centroid2(V,A);
end
end