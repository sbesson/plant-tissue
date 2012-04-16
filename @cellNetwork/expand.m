%% expand
% Expands the cell Network 
%%

%% Syntax   
% N = expand(N,f1,f2)
%
%% Description
% Multiplies the vertices by a fixed number.
%
%% Inputs
% * N - a cellNetwork object
% * f1 - a factor of expansion
% * f2 - an optional parameter
%
%% Outputs
% * N - a new cellNetwork object
%
%% Examples
% >> N = expand(N,2); 
% expand the cellNetwork by a factor 2
%
%% See also 
% * 
%
%% Author 
% Sebastien Besson.
% email address : sbesson@oeb.harvard.edu
% December 2007; Last revision: February 5, 2008

function N = expand(N,f1,f2,delta)

% Check the number of inputs
error(nargchk(2, 4, nargin));

% Test the presence of the optional third argument
if (nargin < 4), delta = 0; end
if (nargin < 3), f2 = f1; end

% Expansion of the vertices
% Radial growth
L=N.v(:,1).^2+ N.v(:,2).^2;
if delta ==1
    N.v(:,1) = N.v(:,1)+N.v(:,1).*(L>=0.98*max(L))*(f1-1);
    N.v(:,2) = N.v(:,2)+N.v(:,2).*(L>=0.98*max(L))*(f2-1);
else
    N.v(:,1) = N.v(:,1)+N.v(:,1)*(f1-1);
    N.v(:,2) = N.v(:,2)+N.v(:,2)*(f2-1);
end

% Isotropic growth
%N.v(:,1) = N.v(:,1)+N.v(:,1).*(f1-1);
%N.v(:,2) = N.v(:,2)+N.v(:,2).*(f2-1);

% Wall stiffening
tau = 10;
N.e(:,4) = N.e(:,4)+(1-N.e(:,4))/tau;

end