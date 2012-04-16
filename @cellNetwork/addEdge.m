%% addEdge
% Method of the cellNetwork class to add an edge
%%

%% Syntax   
% N = addEdge(N,E)
%
%% Description
% Add a new edge to a cellNetwork object.
%
%% Inputs
% * N - a cellNetwork object
% * E - a n by 3 array corresponding to the edges to add
%
%% Outputs
% * N2 - a new cellNetwork object
%
%% Examples
% >> N = addEdge(N,[1 2 0]); 
% add a straight edge between vertices 1 and 2
%
%% See also 
% * removeEdge
%
%% Author 
% Sebastien Besson.
% email address : sbesson@oeb.harvard.edu
% October 2007; Last revision: November 27, 2007

function N = addEdge(N,E)

% Test if the edges can be defined according to the existing vertices
isValid1 = ((max(E(:,1)) <= size(N.v,1)) && (min(E(:,1)) > 0));
isValid2 = ((max(E(:,2)) <= size(N.v,1)) && (min(E(:,2)) > 0));

if (~isValid1 || ~isValid2)
    error('addEdge:error','Edge boundaries are off limits');
end

% Add the edge
switch (size(E,2))
    case 3 
        N.e = [N.e; [E zeros(size(E,1),1)]];
    case 4
        N.e = [N.e; E];
end
% For debugging
%warning('addEdge:operation','Edge %d added',size(N.e,1));
end