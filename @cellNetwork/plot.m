%% plot
% Plot a cellNetwork object.
%%
%% Syntax
% h = plot(N,varargin)
%
%% Description 
%  This
%
%% Inputs
% * N - a cellNetwork object
% * varargin - an optional parameter containing the degree of refinement of
% the plotted edges
%
%% Outputs
% h is a structure containing two arrays of plot handles:
% * h.vertices are the handles
% * h.edges are the handles
%
%% Example
%  >> plot(N,'-r')
%  >> h = plot(N,'-r')
%  >> h = plot(N,'-r','Linewidth',2)
%
%% See also
%
%% Author
% Sebastien Besson
% email address : sbesson@oeb.harvard.edu
% November 2007; Last revision:  November 13, 2008

function h = plot(N,varargin)
h.vertices = plot([],[]);
h.edges = plot([],[]);
h.cells = plot([],[]);

% Detect color specification
if (round(length(varargin)/2)~=(length(varargin)/2))
    % Even number of line parameters - linespec specified 
    linespec=varargin{1};
    % List of colors handled by the linespec argument
    colorlist = {'r';'g';'b';'c';'m';'y';'k';'w'};
    res = regexp(linespec,colorlist,'match','once');
    if ~isempty([res{:}]), color =[res{:}]; end

    plot_arguments = varargin(2:end);
else
    plot_arguments = varargin;
end

for k=1:2:length(plot_arguments) 
    if strcmpi(plot_arguments{k},'color'), color=plot_arguments{k+1}; end
end  
    
hold on
% Plot the vertices
for i = 1:size(N.v,1)
    h.vertices(i) = plot(N.v(i,1), N.v(i,2),varargin{:},'LineStyle','none');
end

% Plot the edges
for i = 1:size(N.e,1)
    if exist('color','var') 
        % If color is specified in the varargin
        newedge = edge(N,i);
        edgecolor = color;
    else
        % If color is not specified in the varargin
        % Retrieve color from the age of the edges
        [newedge,edgecolor] = edge(N,i);
    end
    h.edges(i) = plot(newedge(:,1), newedge(:,2),varargin{:},...
        'Color',edgecolor,'Marker','none');    
end

% Plot the cells
for c = 1:length(N.c)
    celledge = [];
    for j =1:length(N.c{c})
        celledge = [celledge; edge(N,N.c{c}(j))];    
    end
    % Fill the cell with the options specified by varargin
    h.cells(c) = fill(celledge(:,1),celledge(:,2),'k','EdgeColor','none','FaceColor','none');
end

uistack([h.edges],'top');
uistack([h.vertices],'top');

end