%% cell2tri
% Converts a cellNetwork object into a triangle object.
%%

%% Syntax   
% T = cell2tri(N,nc)
%
%% Description
% Create a new triangle object from a cell of a cellNetwork object  
%
%% Inputs
% * N - a cellNetwork object
% * nc - the index of the cell to transform
%
%% Outputs
% * T - a triangle object 
%
%% Examples
% >> T = cell2tri(N,1);
%
%
%% See also 
% * 
%
%% Author 
% Sebastien Besson.
% email address : sbesson@oeb.harvard.edu
% January 2008; Last revision: May 1, 2008

function T = cell2tri(N,nc)

if nargin == 1, nc = 1; end

if (length(N.c{nc}) ~=3)
    error('cell2tr:operation','Dimension of the cell is not 3'); 
end

C = N.c{nc};

A(logical(C>0)) = N.e(C(logical(C>0)),3);
A(logical(C<0)) = -N.e(-C(logical(C<0)),3);

if(cellArea(N,nc)>= 0)
    T = triangle(A(1)+pi/4,A(2)+pi/4);
else
    T = triangle(A(1)-pi/4,A(2)-pi/4);
end
end
