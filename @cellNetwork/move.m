%% move
% Moves the cell Network 
%%

%% Syntax   
% N = move(N)
%
%% Description
% Multiplies the vertices by a fixed number.
%
%% Inputs
% * N - a cellNetwork object
%
%% Outputs
% * N - a new cellNetwork object
%
%% Examples
% >> N = move(N,2); 
%
%% See also 
% * 
%
%% Author 
% Sebastien Besson.
% email address : sbesson@oeb.harvard.edu
% December 2007; Last revision: February 5, 2008

function N = move(N)

% Check the number of inputs
error(nargchk(1, 1, nargin));

% Wall stiffening
%tau = 10;
%N.e(:,4) = N.e(:,4)+(1-N.e(:,4))/tau;

% Calculation of the angles formed between the walls and the x-axis
angles = zeros(size(N.e,1));
for j = 1:size(N.e)
    V=N.v(N.e(j,2),:)-N.v(N.e(j,1),:);
    angles(j) = atan2(V(2),V(1)); %-N.e(j,3);
end

% Radius modification
beta=.2;
for i = 1:size(N.v)
    index1 = find(N.e(:,1) == i);
    index2 = find(N.e(:,2) == i);
    alpha = [angles(index1)-N.e(index1,3); angles(index2)+pi+N.e(index2,3)];
    tensions = [N.e(index1,4); N.e(index2,4)];
    dV(i,:) = beta * [sum(tensions.*cos(alpha))...
         sum(tensions.*sin(alpha))];
end
%dV

% Angle modification
for j = 1:size(N.e,1)
    L0 = N.v(N.e(j,2),:)-N.v(N.e(j,1),:);
    L = L0 + dV(N.e(j,2),:) - dV(N.e(j,1),:);
    area = L0(2)*dV(N.e(j,2),1)-L0(1)*dV(N.e(j,2),2);
    area = area-L(2)*dV(N.e(j,1),1)+L(1)*dV(N.e(j,1),2);
   
    % Test the arc approximation    
    if N.e(j,3)^2 <= 0.00001
       % If the arc approximation is valid
       area = norm(L0)^2*sin(N.e(j,3))/6;
    else
       % General expression for the area
       area = norm(L0)^2/(2*sin(N.e(j,3)))^2 ...
       *(N.e(j,3)-cos(N.e(j,3))*sin(N.e(j,3)));
    end
    %area
    N.e(j,3)= fzero(@(x) 4*area*sin(x)^2-norm(L)^2*(x-cos(x)*sin(x)), N.e(j,3));
end

% Vertex movement
N.v = N.v +dV;

end