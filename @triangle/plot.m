%% plot
% Plot a triangle object in the triangles space.
%%

%% Syntax
% plot(T,delta)
%
%% Description 
%  This function plots the image of the triangle in the space of the
%  triangles.
%
%% Inputs
% * T - a triangle object
% * delta - an optional parameter containing the degree of refinement of
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
% January 2008; Last revision:  April 30, 2008

function h = plot(T,varargin)

% Check the number of inputs

x = abs([T.a1]);
y = abs([T.a2]);
hold on

% Scott's space of triangle
%alpha=atan(.5);
%X = (pi/2-y)*sin(alpha)-(x-pi/4)*cos(alpha);
%Y = (pi/2-y)*cos(alpha)+(x-pi/4)*sin(alpha);

% Simpler space of triangles
X = y/cos(pi/6)+x*tan(pi/6);
Y = x;

% Plots the vertices
if (nargin > 1)
    h = plot(X,Y,varargin{:});
else
    h = plot(X,Y);
end

axis tight equal off;
end