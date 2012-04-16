%% rotate
% Rotates the cell Network 
%%

%% Syntax   
% N = rotate(N,C,angle)
%
%% Description
% Applies a rotation to the cellNetwork object.
%
%% Inputs
% * N - a cellNetwork object
% * C - the center of rotation
% * angle - the rotation angle 
%
%% Outputs
% * N - a new cellNetwork object
%
%% Examples
% >> N = rotate(N,[0 0],pi/3); 
% rotate the cellNetwork with the vector [1,2]
%
%% See also 
% * 
%
%% Author 
% Sebastien Besson.
% email address : sbesson@oeb.harvard.edu
% March 2008; Last revision: March 29, 2008

function N = rotate(N,C,angle)

% Check the number of inputs
error(nargchk(3, 3, nargin));

% Check the dimension of the translation vector
if (length(C) ~= 2)
    error('rotate:dimensions',...
        'Dimensions of the translation vector must be equal to 2.') 
end
    
% Rotates the vertices
R = [cos(angle) -sin(angle); sin(angle) cos(angle)];
N.v = (N.v - repmat(C,size(N.v,1),1))*R + repmat(C,size(N.v,1),1);
    
end