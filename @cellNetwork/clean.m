%% clean
% Cleans a cellNetwork object
%%

%% Syntax   
% N = clean(N)
%
%% Description
% Loops over the defined cells and remove unused vertices and edges.
%
%% Inputs
% * N - a cellNetwork object
%
%% Outputs
% * N - a new cellNetwork object
%
%% Examples
% >> N = clean(N); 
% cleans the cellNetwork object
%
%% See also 
% * removeEdge
% * removeVertex
%
%% Author 
% Sebastien Besson.
% email address : sbesson@oeb.harvard.edu
% October 2007; Last revision: November 27, 2007

function N = clean(N)

% Remove cells that are not well defined
for i=1:length(N.c)
    Vstart = zeros(length(N.c{i}),2);
    Vend = zeros(length(N.c{i}),2);
    Vstart(N.c{i}>0,:) = N.v(N.e(N.c{i}(N.c{i}>0),1),:);
    Vstart(N.c{i}<0,:) = N.v(N.e(-N.c{i}(N.c{i}<0),2),:);
    Vend(N.c{i}>0,:) = N.v(N.e(N.c{i}(N.c{i}>0),2),:);
    Vend(N.c{i}<0,:) = N.v(N.e(-N.c{i}(N.c{i}<0),1),:);
    if any(any(Vstart-circshift(Vend,[1 0])))
        removeCell(N,i);
    end
end

% List existing edges
ne = 1:size(N.e,1);

% Remove edges belonging to defined cells
ne(unique(abs([N.c{:}]))) = [];

% Remove unused edges
if ~isempty(ne), N = removeEdge(N,ne); end;


% List existing vertices
nv = 1:size(N.v,1);

% Remove vertices belonging to defined cells
nv(unique(N.e(:,1:2))) = [];

% Remove unused vertices
if ~isempty(nv), N = removeVertex(N,nv); end;

end