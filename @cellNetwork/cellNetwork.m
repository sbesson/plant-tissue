%% cellNetwork
% cellNetwork class constuctor
%%

%% Syntax   
% N = cellNetwork(V,E,C)
%
%% Description
% Create a new object of the type cellNetwork.
%
%% Inputs
% * V - a n by 2 matrix containing the coordinates of the vertices
% * E - a m by 3 matric containing the edge information
% * C - a l by 1 array of cells containing row vectors defining the cells
% from their edges
%
%% Outputs
% * N - a cellNetwork object 
%
%% Examples
% >> N =  cellNetwork();
% creates an empty cellNetwork object.
% >> N2 =  cellNetwork(N);
% copies a cellNetwork object.
% >> N =  cellNetwork(V,E,C);
% creates a cellNetwork object with V as vertices, E as edges and C cells.
%
%% See also 
% * 
%
%% Author 
% Sebastien Besson.
% email address : sbesson@oeb.harvard.edu
% October 2007; Last revision: November 26, 2007

function N = cellNetwork(varargin)

switch nargin
    case 0
        % if no input arguments, create a default object
        N.v = [];
        N.e= [];
        N.c = {};
        N = class(N,'cellNetwork');
    case 1
        % if single argument of class cellNetwork, return it 
        if (isa(varargin{1},'cellNetwork')) 
            N = varargin{1};
        else
            error('Wrong argument type')  
        end 
    case 3
        % create object using specified values
        N.v = varargin{1};
        if (size(varargin{2},2) == 3)
            N.e = [varargin{2} ones(size(varargin{2},1),1)];
        else
            N.e = varargin{2};
        end
        N.c = varargin{3};
        N = class(N,'cellNetwork');
    otherwise
         error('Wrong number of arguments') 
end
check(N)
    
end