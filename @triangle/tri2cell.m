%% tri2cell
% Converts a triangle object into a cellNetwork object.
%%

%% Syntax   
% N = tri2cell(T)
%
%% Description
% Create a new cellNetwork  
%
%% Inputs
% * T - a triangle object
%
%% Outputs
% * N - a cellNetwork object 
%
%% Examples
% >> N = tri2cell(triangle(pi/3,pi/2));
%
%
%% See also 
% * cell2tri 
%
%% Author 
% Sebastien Besson.
% email address : sbesson@oeb.harvard.edu
% January 2008; Last revision: April 10, 2008

function N = tri2cell(T,angle)

if (nargin == 1), angle = pi/2; end;

x = tan(T.a2)/(tan(T.a1)+tan(T.a2));
y = tan(T.a1)*tan(T.a2)/(tan(T.a1)+tan(T.a2));

%Vertices
V = [0 0; 
    1 0;
    x y];
theta = atan2(T.V2(2)-T.V1(2),T.V2(1)-T.V1(1));
V = norm(T.V2-T.V1)*V*[cos(theta) sin(theta);-sin(theta) cos(theta)]...
    +repmat(T.V1,3,1);
% Edges
E = [2 3 sign(T.a1)*(abs(T.a1)-pi/4+(angle-pi/2)) 1;
     3 1 sign(T.a1)*(abs(T.a2)-pi/4) 2;
     1 2 sign(T.a1)*(3*pi/4-abs(T.a2)-abs(T.a1)) 3;];
% Cells
C = {[1 2 3]};
N = cellNetwork(V,E,C);
end