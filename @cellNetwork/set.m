%% set
% cellNetwork class set method
%%

%% Syntax   
% val = set(N,PropName)
%
%% Description
% Sets the property of a cellNetwork object.
%
%% Inputs
% * N - a cellNetwork object 
% * PropName - a string containing the propertu to change
%
%% Outputs
%
%% Examples
% >> N = set(N,'Vertices',V);
% >> N = set(N,'Edges',E);
% >> N = set(N,'Cells',C);
%
%% See also 
% * get
%
%% Author 
% Sebastien Besson.
% email address : sbesson@oeb.harvard.edu
% October 2007; Last revision: November 26, 2007

function N = set(N,PropName,val)

switch PropName
    case 'Vertices'
        N.v = val;
    case 'Edges'
        N.e = val;
    case 'Cells'       
        N.c = val;
    otherwise
        error([propName ,'Is not a valid stock property']);
end
end