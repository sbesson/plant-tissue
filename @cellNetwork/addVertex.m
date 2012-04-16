%% addVertex
% Method of the cellNetwork class to add a Vertex
%%

%% Syntax   
% N = addVertex(N,V)
%
%% Description
% Add a new vertex to a cellNetwork object.
%
%% Inputs
% * N - a cellNetwork object
% * V - a n by
%
%% Outputs
% * N - a new cellNetwork object
%
%% Examples
% >> N = addVertex(N,[0 0]); 
% add the origin to the N network
%
%% See also 
% * addEdge
%
%% Author 
% Sebastien Besson.
% email address : sbesson@oeb.harvard.edu
% October 2007; Last revision: November 26, 2007

function N = addVertex(N,V)

N.v = [N.v; V];

end