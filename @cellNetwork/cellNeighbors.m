%% cellNeighbors
% Find neighbors of a series of cells in a cellNetwork object.
%%

%% Syntax
% NC = cellNeighbors(N,n)

%% Description 
%  This
%
%% Inputs
% * N - a cellNetwork object
% * n - an array of cell indices
%
%% Outputs
% * NC - a cell array containing the neighbor cells
%
%% Example
%
%% See also
% * neighbors
%
%% Author
% Sebastien Besson
% email address : sbesson@oeb.harvard.edu
% April 2008; Last revision:  April 13, 2009

function ncells = cellNeighbors(N,nc)

% Check the number of inputs
error(nargchk(1, 2, nargin));

% If no optional argument, loops over the cells
if (nargin == 1), nc = 1:size(N.c,1); end;

% Initialize the output array
for i=1:length(nc)
    ncells{i} = unique(neighbors(N,abs(N.c{nc(i)})));
end

end