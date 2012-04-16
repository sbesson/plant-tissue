%% check
% Check a cellNetwork object
%%

%% Syntax   
% check(N)
%
%% Description
% Loops over the defined cells and checks unused vertices and edges as well as badly defined cells.
%
%% Inputs
% * N - a cellNetwork object
%
%% Outputs
%
%% Examples
% >> check(N); 
% check validity of the cellNetwork object
%
%% See also 
% * clean
%
%% Author 
% Sebastien Besson.
% email address : sbesson@oeb.harvard.edu
% October 2009; Last revision: October 23, 2009

function  check(N)

% Remove cells that are not well defined
for i=1:length(N.c)
    Vstart = zeros(length(N.c{i}),2); % List of start vertices
    Vend = zeros(length(N.c{i}),2); % List of end vertices
    Vstart(N.c{i}>0,:) = N.v(N.e(N.c{i}(N.c{i}>0),1),:);
    Vstart(N.c{i}<0,:) = N.v(N.e(-N.c{i}(N.c{i}<0),2),:);
    Vend(N.c{i}>0,:) = N.v(N.e(N.c{i}(N.c{i}>0),2),:);
    Vend(N.c{i}<0,:) = N.v(N.e(-N.c{i}(N.c{i}<0),1),:);
    if any(any(Vstart-circshift(Vend,[1 0])))
        error('check:error','Bad definition of cell %g',i);
    end
end

% List existing edges
ne = 1:size(N.e,1);

% Remove edges belonging to defined cells
ne(unique(abs([N.c{:}]))) = [];

% Remove unused edges
if ~isempty(ne), warning('check:warning',['Unused edges:' sprintf('%d ',ne)]); end;

% List existing vertices
nv = 1:size(N.v,1);

% Remove vertices belonging to defined cells
if ~isempty(N.e), nv(unique(N.e(:,1:2))) = []; end

% Remove unused vertices
if ~isempty(nv), warning('check:warning',['Unused vertices:' sprintf('%d ',nv)]); end;

end