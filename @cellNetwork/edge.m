%% edge
% Get the data corresponding to the edge.
%%

%% Syntax
% edge(N,n)
%
%% Description 
%  This
%
%% Inputs
% * N - a network structure
% * n - an optional parameter containing the degree of refinement of
% the plotted edges
%
%% Outputs
% none
%
%% Example
%
%% See also
%
%% Author
% Sebastien Besson
% email address : sbesson@oeb.harvard.edu
% April 2008; Last revision:  May 24, 2008

function [edgedata,color] = edge(N,n,s)

% Default value for the discretization of the edges
if nargin<3
    s = 0:0.02:1;
end

% Test if edge is positive of negatively oriented
if n>0
    edgedata = findP([N.v(N.e(n,1),:);N.v(N.e(n,2),:)],N.e(n,3),s);
    color = coloredge(N.e(n,4));
else
    edgedata = findP([N.v(N.e(-n,2),:);N.v(N.e(-n,1),:)],-N.e(-n,3),s);
    color = coloredge(N.e(-n,4));
end
end