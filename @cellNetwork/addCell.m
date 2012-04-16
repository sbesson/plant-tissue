%% addCell
% Method of the cellNetwork class to add a cell
%%

%% Syntax   
% N = addCell(N,C)
%
%% Description
% Add a new cell to a cellNetwork object.
%
%% Inputs
% * N - a cellNetwork object
% * C - a new cell
%
%% Outputs
% * N2 - a new cellNetwork object
%
%% Examples
% >> N = addCell(N,[1 2 4]); 
% add a cell defined by the vertices 1, 2 and 4
%
%% See also 
% * removeCell
%
%% Author 
% Sebastien Besson.
% email address : sbesson@oeb.harvard.edu
% October 2007; Last revision: May 14, 2009

function N = addCell(N,C)

% Add the cell
N.c{end+1} = C;

% For debugging
warning('addCell:operation','Cell %d added',size(C,1));
end