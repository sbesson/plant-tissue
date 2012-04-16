%% cellBoundary
% Find the intersection of a half line with a cell boundary.
%%

%% Syntax   
% B = cellBoundary(N,n,angles);
%
%% Description 
% This program computes the intersection of the line starting at the
% centroid of the cell and forming an angle theta with the horizontal line
% with the boundary of the cell. 
%
%% Inputs
% * N - a celNetwork object
% * theta - the angle formed by the line with the horizontal axis in the
% range [0, 2*pi]
% 
%% Outputs
% * B - the boundary network
%
%% Examples
% >> cellBoundary(N,0:pi/100:pi/4)
%
%% TODO
%
%% See also
% * findP
% * cellcentre
%
%% Author    
% Sebastien Besson
% email address: sbesson@oeb.harvard.edu
% April 2008; Last revision: May 28 2007

function B = cellBoundary(N,n,angles)

% Check number of inputs
error(nargchk(3, 3, nargin));

% If no optional argument, loops over the cells
if (nargin == 1), n = 1:length(N.c); end;

B = zeros(length(angles),2,length(n));
for c = 1:length(n)    
    %Extracts the cell data
    C = N.c{n(c)};
    
    % Extracts the vertices and angles
    V = zeros(length(C),2);
    A = zeros(length(C),1);
    
    index1 = logical(C>0);
    V(index1,:) = N.v(N.e(C(index1),1),:);
    A(index1) = N.e(C(index1),3); 
    
    index2 = logical(C<0);
    V(index2,:) = N.v(N.e(-C(index2),2),:);
    A(index2) = - N.e(-C(index2),3); 


     B(:,:,n(c))= boundary(V,A,angles);
end

end