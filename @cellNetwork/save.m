%% save
% Save a cellNetwork object as an SVG file.
%%
%% Syntax
% save(N,filename,varargin)
%
%% Description 
% This function saves the cellNetwork as a collection of edges using the
% SVG format.
%
%% Inputs
% * N - a cellNetwork structure
% * filename - the filename under which the cellNetwork will be saved
% * varargin - an optional parameter
%
%% Outputs
%
%% Example
%  >> save(N,'object.svg')
%
%% See also
%
%% Author
% Sebastien Besson
% email address : sbesson@oeb.harvard.edu
% November 2007; Last revision:  October 21, 2008


function save(N,filename,varargin)

fid=fopen(filename, 'wb');

if fid <= 0
  error('failed to open file');
end

fprintf(fid, '%s\n', '<?xml version="1.0" standalone="no"?>');
fprintf(fid, '%s\n', '<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" ');
fprintf(fid, '%s\n', '"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">');
fprintf(fid, '%s\n', '<svg preserveAspectRatio="xMinYMin meet" width="100%" height="100%" version="1.1" xmlns="http://www.w3.org/2000/svg">');

% Set default values
stroke='black';
stroke_width=1;

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

% parameter handling
for argindex=1:2:length(plot_arguments)
  arg=varargin{argindex};
  param=varargin{argindex+1};
  switch(arg)
   case 'Color'
    stroke=param;
   case 'Linewidth'
    stroke_width=param;
   otherwise
    error(['invalid parameter at position ', mat2str(argindex)]);
  end % switch
end % for argindex

% Loop over edges
for i =1:size(N.e,1)
    fprintf(fid, '<path d="');
    % initial point
    fprintf(fid, 'M %f %f ', 100*N.v(N.e(i,1),:));
    if (abs(N.e(i,3))<=0.001)
        % If edge is a line
        fprintf(fid, 'L %f %f ', 100*N.v(N.e(i,2),:));
    else
        % If edge is a curve
        R = 100*norm(N.v(N.e(i,2),:)- N.v(N.e(i,1),:))/(2*sin(N.e(i,3)));
        ispositive = (N.e(i,3)>0);
        fprintf(fid, 'A %f %f 0 0 %f %f %f ', R,R,ispositive,100*N.v(N.e(i,2),:));
    end
    if ~exist('color','var'),
        ecolor = coloredge(N.e(i,4));
        stroke = sprintf('rgb(%.0f,%.0f,%.0f)',ecolor*255);
    end
    fprintf(fid, '%s', '" style="fill:none;');
    fprintf(fid, 'stroke:%s;stroke-width:%i"/>', stroke, stroke_width);
    fprintf(fid, '\n');
end
fprintf(fid, '</svg>\n');
fclose(fid);
end