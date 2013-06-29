%% sanityCheck
% Check a cellNetwork object
%%

%% Syntax   
% sanityCheck(N)
%
%% Description
% Loops over the defined cells and sanityChecks unused vertices and edges as well as badly defined cells.
%
%% Inputs
% * N - a cellNetwork object
%
%% Outputs
%
%% Examples
% >> sanityCheck(N); 
% sanityCheck validity of the cellNetwork object
%
%% See also 
% * clean
%
%% Author 
% Sebastien Besson.
% email address : sbesson@oeb.harvard.edu
% October 2009; Last revision: October 23, 2009

function  sanityCheck(N)

% Remove cells that are not well defined
nCells = length(N.c);
for i = 1 : nCells
    Vstart = zeros(length(N.c{i}),2); % List of start vertices
    Vend = zeros(length(N.c{i}),2); % List of end vertices
    Vstart(N.c{i}>0,:) = N.v(N.e(N.c{i}(N.c{i}>0),1),:);
    Vstart(N.c{i}<0,:) = N.v(N.e(-N.c{i}(N.c{i}<0),2),:);
    Vend(N.c{i}>0,:) = N.v(N.e(N.c{i}(N.c{i}>0),2),:);
    Vend(N.c{i}<0,:) = N.v(N.e(-N.c{i}(N.c{i}<0),1),:);
    if any(any(Vstart-circshift(Vend,[1 0])))
        error('sanityCheck:error','Bad definition of cell %g',i);
    end
end

% List existing edges
ne = 1:size(N.e,1);

% Remove edges belonging to defined cells
for i = 1 : nCells, ne(abs(N.c{i})) = NaN; end
ne(isnan(ne))=[];

% Remove unused edges
if ~isempty(ne), warning('sanityCheck:warning',['Unused edges:' sprintf(' %d',ne)]); end;

% List existing vertices
nv = 1:size(N.v,1);

% Remove vertices belonging to defined cells
if ~isempty(N.e), nv(unique(N.e(:,1:2))) = []; end

% Remove unused vertices
if ~isempty(nv), warning('sanityCheck:warning',['Unused vertices:' sprintf(' %d',nv)]); end;

end