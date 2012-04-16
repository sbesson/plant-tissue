%% replot
% Plot a cellNetwork object.
%%
%% Syntax
% h = replot(N,h,varargin)
%
%% Description 
%  This
%
%% Inputs
% * N - a cellNetwork object
% * h - a handle structure to the plot objects
% * varargin - an optional parameter containing the degree of refinement of
% the plotted edges
%
%% Outputs
% none
%
%% Example
%  >> replot(N,h)
%  >> h = replot(N,h)
%  >> h = replot(N,h,'-r')
%  >> h = replot(N,h,'-r','Linewidth',2)
%
%% See also
%
%% Author
% Sebastien Besson
% email address : sbesson@oeb.harvard.edu
% November 2007; Last revision:  November 13, 2008

function h = replot(N,h,varargin)

if ~isempty(h.vertices)
    parentobject = get(h.vertices(1),'Parent');
else
    parentobject = gca;
end
% Detect color specification
if (round(length(varargin)/2)~=(length(varargin)/2))
    % Even number of line parameters - linespec specified
    linespec=varargin{1};
    % List of colors handled by the linespec argument
    colorlist = {'r';'g';'b';'c';'m';'y';'k';'w'};
    res = regexp(linespec,colorlist,'match','once');
    if ~isempty([res{:}]), color =[res{:}]; end
    
    plot_arguments = varargin(2:end);
    % List of linestyles handled by the linespec argument
    linelist = {'-';'--';':';'_.'};
    res = regexp(linespec,linelist,'match','once');
    if ~isempty([res{:}])
        plot_arguments = {plot_arguments{:} 'LineStyle' [res{:}]}; 
    end

    % List of markerstyles handled by the linespec argument
    markerlist = {'+';'o';'*';'\.';'x';'s';'d';'^';'v';...
        '<';'>';'p';'h'};
    res = regexp(linespec,markerlist,'match','once');
    if ~isempty([res{:}])
        plot_arguments = {plot_arguments{:} 'Marker' [res{:}]};
    end
else
    plot_arguments = varargin;
end

k=1;
while k<length(plot_arguments) 
    if strcmpi(plot_arguments{k},'color'), color=plot_arguments{k+1}; end
    k= k+2;
end  

% Test if number of vertices has decreased
if length(h.vertices)>size(N.v,1)
    % Delete addditional vertices
    delete(h.vertices(size(N.v,1)+1:length(h.vertices)));
    h.vertices(size(N.v,1)+1:length(h.vertices))=[];
end
% Test if the number edges has decreased
if length(h.edges)>size(N.e,1)
    % Delete additional edges
    delete(h.edges(size(N.e,1)+1:length(h.edges)));
    h.edges(size(N.e,1)+1:length(h.edges))=[];
end

% Test if the number edges has decreased
if length(h.cells)>length(N.c)
    % Delete additional edges
    delete(h.cells(length(N.c)+1:length(h.cells)));
    h.cells(length(N.c)+1:length(h.cells))=[];
end

hold on
% Replot vertices 
for i = 1:size(N.v,1)
    if i<=length(h.vertices)
        % If graphic vertex already exist
        set(h.vertices(i),'XData',N.v(i,1),'YData',N.v(i,2),...
            plot_arguments{:},'LineStyle','none');
    else
        % Add new graphic vertex
        h.vertices(i) = plot(N.v(i,1),N.v(i,2),...
            varargin{:},'LineStyle','none','Parent',parentobject);
    end
end

% Replot  edges 
for i = 1:size(N.e,1)
    if exist('color','var') 
        % If color is specified in the varargin
        newedge = edge(N,i);
        edgecolor=color;
    else
        % If color is not specified in the varargin
        % Retrieve color from the age of the edges
        [newedge,edgecolor] = edge(N,i);
    end
    
    if i<=length(h.edges)
        % If graphic edge already exist
        set(h.edges(i),'XData',newedge(:,1),'YData',newedge(:,2),...
            plot_arguments{:},'Color',edgecolor,'Marker','none');
    else
        % Add new graphic edge
        h.edges(i) = plot(newedge(:,1), newedge(:,2),varargin{:},...
            'Color',edgecolor,'Marker','none','Parent',parentobject);
    end
end

% Replot  cells 
for c = 1:length(N.c)

    celledge = [];
    for j =1:length(N.c{c})
        celledge = [celledge; edge(N,N.c{c}(j))];    
    end
    
    if c<=length(h.cells)
        set(h.cells(c),'XData',celledge(:,1),'YData',celledge(:,2));
    else
        % Fill the cell with the options specified by varargin
        h.cells(c) = fill(celledge(:,1),celledge(:,2),'k','EdgeColor','none','FaceColor','none');
    end
    hold on;
end


uistack([h.vertices],'top')
end