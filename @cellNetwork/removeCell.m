%% removeCell
% Method of the cellNetwork class to remove a cell
%%

%% Syntax   
% N = removeCell(N,nc)
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
% >> N = removeCell(N,1); 
% remove the first cell of the network
%
%% See also 
% * addCell
%
%% Author 
% Sebastien Besson.
% email address : sbesson@oeb.harvard.edu
% October 2007; Last revision: December 6, 2007

function N = removeCell(N,nc)

% Test if the cells exists
if (min(nc) <= 0)
    error('removeCell:arraySize','Cell %d does not exist',min(nc));
end

if (max(nc) > length(N.c))
    error('removeCell:arraySize','Cell %d does not exist',max(nc));
end

% Removes the cells
N.c(nc) =[];

% For debugging
for i = 1:size(nc,2)
    warning('removeCell:operation','Cell %d removed',nc(i));
end

end