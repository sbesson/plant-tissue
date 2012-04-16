%% get
% cellNetwork class get method
%%

%% Syntax   
% val = get(N,PropName)
%
%% Description
% Return the property of a cellNetwork object.
%
%% Inputs
% * N - a cellNetwork object 
% * PropName - a string containing the Information
%
%% Outputs
% * N - a cellNetwork object 
%
%% Examples
% >> V = get(N,'Vertices')
% >> E = get(N,'Edges')
% >> C = get(N,'Cells')
%
%% See also 
% * set
%
%% Author 
% Sebastien Besson.
% email address : sbesson@oeb.harvard.edu
% October 2007; Last revision: November 26, 2007

function val = get(N,PropName)

switch PropName
    case 'Vertices'
        val = N.v;
    case 'Edges'
        val = N.e;
    case 'Cells'       
        val = N.c;
    otherwise
        error([PropName,' is not a valid property']);
end
end