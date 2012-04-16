%% edgelength
% Get the length corresponding to the edge.
%%

%% Syntax
% edgelength(N,n)
%
%% Description 
%  This
%
%% Inputs
% * N - a network structure
% * n - a vector containing the list of edges
%
%% Outputs
% none
%
%% Example
%
%% See also
%
%% Author
% Sebastien Besson
% email address : sbesson@oeb.harvard.edu
% May 2008; Last revision:  May 8, 2009

function L = edgelength(N,n)

% Check the number of inputs
error(nargchk(1, 2, nargin));

% If no optional argument, loops over the cells
if (nargin == 1), n = 1:size(N.e,1); end;

V1 = N.v(N.e(abs(n),1),:);
V2 = N.v(N.e(abs(n),2),:);

dV = sqrt(sum((V2-V1).^2,2));
A=N.e(abs(n),3);
smallA = A.^2<=0.00001;
largeA = A.^2>0.00001;
L(smallA) = dV(smallA).*(1+1/6*A(smallA).^2);
L(largeA) = dV(largeA).*A(largeA)./sin(A(largeA));

end