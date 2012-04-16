%% translate
% Translates the cell Network 
%%

%% Syntax   
% N = translate(N,T)
%
%% Description
% Multiplies the vertices by a fixed number.
%
%% Inputs
% * N - a cellNetwork object
% * T - a translation vector
%
%% Outputs
% * N - a new cellNetwork object
%
%% Examples
% >> N = translate(N,[1 2]); 
% translate the cellNetwork with the vector [1,2]
%
%% See also 
% * 
%
%% Author 
% Sebastien Besson.
% email address : sbesson@oeb.harvard.edu
% January 2008; Last revision: March 29, 2008

function N = translate(N,T)

% Check the number of inputs
error(nargchk(2, 2, nargin));

% Check the dimension of the translation vector
if (length(T) ~= 2)
    error('translate:dimensions',...
        'Dimensions of the translation vector must be equal to 2.') 
end
    
% Translates the vertices
N.v = N.v + repmat(T,size(N.v,1),1);
    
end