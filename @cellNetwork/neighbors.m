%% neighbors
% Find neighbors cells with respect to an edge a cellNetwork object.
%%

%% Syntax
% NC = neighbors(N,n)

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
%
%% Author
% Sebastien Besson
% email address : sbesson@oeb.harvard.edu
% June 2008; Last revision:  December 3, 2008

function ncells = neighbors(N,ne)

% Check the number of inputs
error(nargchk(1, 2, nargin));

% If no optional argument, loops over the edges
if (nargin == 1), ne = 1:size(N.e,1); end;

% Initialize the output array
ncells = [];

for c=1:length(N.c)
    if ~isempty(intersect(abs(ne),abs(N.c{c}))), ncells = [ncells c]; end;        
end
%ncells=zeros(length(ne),2)
%for i=1:length(ne)
%    i
%     %cellfun(@(x) display(x),N.c);
%     [x,y,z] = cellfun(@(x) intersect(ne(i),abs(x)),N.c,'UniformOutput',false);
%     z
%     [z{:}]
%     if ~isempty(length([z{:}]))
%        ncells(i,[z{:}]) = [z{:}]
%     end
%end

end