%% merge
% Merge two cells of a network.
%%

%% Syntax
%  N = merge(N,n1,n2);
%
%% Description
% Merge two given cells of a cellNetwork object
%
%% Inputs
% * N - a network structure
% * n1,n2 - the indices of two cells
%
%% Outputs
% * N - the new network structure
%% See also
% * % divide
%
%% TODO
% to be verified
%
%% Author
% Sebastien Besson
% email address: sbesson@oeb.harvard.edu
% December 2008; Last revision: December 3, 2008

function N = merge(N,n)

% Check thenumber of inputs
error(nargchk(2, 2, nargin));

% Find the common edges of both cells and their positions
[e i1 i2] = intersect(abs(N.c{n(1)}),abs(N.c{n(2)}));
ne = length(e);
nc1= length(N.c{n(1)});
nc2= length(N.c{n(2)});

% Shift the indices so that the cell is of type
% [1:ne ne+1:end] 
M = sort(mod(repmat(i1,nc1,1)-repmat((1:nc1)',1,length(i1)),nc1)+1,2);
[junk d] = intersect(M,1:ne,'rows');
cell1 = circshift(N.c{n(1)},[1 -d+1]);

% Shift the indices so that the cell is of type
% [1:ne ne+1:end] 
M = sort(mod(repmat(i2,nc2,1)-repmat((1:nc2)',1,length(i2)),nc2)+1,2);
[junk d] = intersect(M,1:ne,'rows');
cell2 = circshift(N.c{n(2)},[1 -d+1]);

% Define oriented parts of cells to stick together
newcell1= cell1(ne+1:end);
if all(cell1(1:ne) == cell2(1:ne))
    newcell2= -cell2(end:-1:ne+1);
elseif all(cell1(1:ne) == -cell2(ne:-1:1))
    newcell2= cell2(ne+1:end);
else
    error('cellNetwork::merge');
end

% Add the new cell
N=addCell(N,[newcell1 newcell2]);

% Remove the old cells
N=removeCell(N,n);

% Remove the edges
N=clean(N);

end