%% extractCell
% Method of the cellNetwork class to extract a series of cells
%%

%% Syntax   
% N = extractCell(N,nc)
%
%% Description
% Remove a cell from a cellNetwork object.
%
%% Inputs
% * N - a cellNetwork object
% * nc - a row vector containing the indices of the cells to be removed
%
%% Outputs
% * N - a new cellNetwork object
%
%% Examples
% >> N = extractCell(N,1); 
% extract the first cell of the network
% >> N = extractCell(N,[3 6]); 
% extract the cells 3 and 6 of the network
%
%% See also 
% * addCell
%
%% Author 
% Sebastien Besson.
% email address : sbesson@oeb.harvard.edu
% October 2007; Last revision: December 6, 2007

function N = extractCell(N,nc)

% Test if the cells exists
if (min(nc) <= 0)
    error('removeCell:arraySize','Cell %d does not exist',min(nc));
end

if (max(nc) > length(N.c))
    error('removeCell:arraySize','Cell %d does not exist',max(nc));
end

% Removes the cells
ic = 1:length(N.c);
ic(nc) = [];
N.c(ic) =[];

N = clean(N);
end