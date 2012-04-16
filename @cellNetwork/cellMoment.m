%% cellMoment
% Computes the moments of area of a cellNetwork object.
%%

%% Syntax
% [S,C,I] = cellMoment(N,nc)
% 
%% Description
% Compute the area, first and second moments of area  of a given array of cells 
%
%% Inputs
% * N - a cellNetwork object
% * n - an array of cell indices
%
%% Outputs
% * S - the area.
% * C - the first moment of area.
% * I - the second moment of area.
%
%
%% Examples
% >> S =  cellMoment(N,1);
% >> [S, C, I] =  cellMoment(N);
%
%% See also
% * cellArea
%
%% Author
% Sébastien Besson
% email address : jdumais@oeb.harvard.edu
% November 2007; Last revision: May 15, 2009

function [S,C,I] = cellMoment(N,n)

% Check the number of inputs
error(nargchk(1, 2, nargin));

% If no optional argument, loops over the cells
if (nargin == 1), n = 1:length(N.c); end;

% Initialize the output array
S = zeros(length(n),1);
C = zeros(length(n),2);
I = zeros(length(n),2,2);

for c = 1:size(n,2)
    nc = n(c);
    
    %Extract the cell data
    edges = N.c{nc};
    
    % Extract the vertices and angles
    V = zeros(length(edges),2);
    A = zeros(length(edges),1);
    
    index1 = logical(edges>0);
    V(index1,:) = N.v(N.e(edges(index1),1),:);
    A(index1) = N.e(edges(index1),3); 
    
    index2 = logical(edges<0);
    V(index2,:) = N.v(N.e(-edges(index2),2),:);
    A(index2) = - N.e(-edges(index2),3); 

    % Compute the moments of area of the cell
    [S(c),C(c,:),I(c,:,:)] = moment(V,A);
end
end