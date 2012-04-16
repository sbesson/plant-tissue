%% dimCell
% Method of the cellNetwork class to get the dimensions of the cells
%%

%% Syntax   
% d = dimCell(N)
%
%% Description
% 
%
%% Inputs
% * N - a cellNetwork object
%
%% Outputs
% * d - a row vector
%
%% Examples
% >> d = dimCell(N,n); 
% returns the dimension of the cells of the cell network N
%
%% See also 
% * addCell
% * removeCell
%
%% Author 
% Sebastien Besson.
% email address : sbesson@oeb.harvard.edu
% October 2007; Last revision: June 24, 2008

function d = dimCell(N,n)

if (isempty(N.c)), error('dimCell:operation','No cell'); end;

% If no optional argument, loop over all cells
if (nargin == 1), n = 1:length(N.c); end;

% Apply length operation to all cells
d = cellfun(@length,N.c(n));

end