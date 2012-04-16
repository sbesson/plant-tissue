%% findLength2
% Compute the length of the dividing wall of a cell.
%%

%% Syntax   
% [length,P,Q,angle] = findLength(s,V,A,i,j,area)
%
%% Description 
% Given the boundary of a cell, this function computes the length of the
% wall dividing the cell between two given positions and satisfying the
% area constraint.
%
%% Inputs
% * s is a 2 element-vector which containing the position of the points
% * V is a n-by-2 matrix with the coordinates of the vertices. 
% * A is the curvature angle of the edges.  A positive (negative) angle
% means that the wall bows outward (inward).
% * i and j - the indexes of the starting and ending edges of the new wall.
% * area - the surface constraint for the new cell.
%
%% Outputs
% * length- a variable containing the length of the new wall.
% * P - the coordinates of the first point.
% * Q - the coordinates of the second point.
% * angle - the angle of the new wall.
%
%% Examples
% >> L =  findLength([0.3 0.7],[0 0;0 1; 1 1; 1 0], [0 0 0 0],2,3,pi/8);
%
%% See also
% * findP
% * cellArea
%
%% Author
% Sebastien Besson
% email address : sbesson@oeb.harvard.edu
% October 2007; Last revision: November 12 2007

function [length,P,Q,angle] = findLength2(s1,s2,V,A,i,j,area)

if(i ==j)
    s =sort([s1 s2]);
    s1 = s(1);
    s2 = s(2);
end
    
% Extract s1 and s2 and norms these values between 0 and 1
%s1 = min(max(s1,0),1);
%s2 = min(max(s2,0),1);

% Find the coordinates of the 2 new points on the edges
P = findP(V(i:i+1,:),A(i),s1);
Q = findP(V(j:j+1,:),A(j),s2);

% Define the daughter cell
if (j>i)
    % The daughter cell is[P V(i+1) ... V(j) Q]
    newV = [P; V(i+1:j,:); Q];
    newA = [A(i)*(1-s1); A(i+1:j-1); A(j)*s2; 0];
elseif (i>j)
    % The daughter cell is[P V(i+1)...V(size(V)) V(1)...V(j) Q]
    newV = [P; V(i+1:end-1,:); V(1:j,:); Q];
    newA = [A(i)*(1-s1); A(i+1:end); A(1:j-1); A(j)*s2; 0];
else
    % The half cell is [P Q]
    newV = [P; Q];
    newA = [A(i)*(s2-s1); 0];
end

% Determines the area of the (newV,newA)
cArea=0;
newV = [newV;newV(1,:)];
for nv = 1:size(newV,1)-1
    L = norm(newV(nv+1,:)-newV(nv,:));
    % Test the arc approximation
    if newA(nv)^2 <= 0.00001
       % If the arc approximation is valid
       subArea = L^2*sin(newA(nv))/6;
    else
       % General expression for the area
       subArea = L^2*(newA(nv)-cos(newA(nv))*sin(newA(nv)))/(2*sin(newA(nv)))^2;
    end

    cArea = cArea + (newV(nv,1)*newV(nv+1,2)-newV(nv+1,1)*newV(nv,2))/2 + subArea;
end

subarea = area - cArea;
L= norm(Q-P);

% To satisty the area condition, the angle of the cord [QP] must satisty
% the equation

% Find the angle of the wall satisfying the area condition  
angle= fzero(@(x) 4*subarea*sin(x)^2-L^2*(x-cos(x)*sin(x)), 1);
%if (angle <= 0.000001), angle = fzero(@(x) 6*subarea-L^2*sin(x),angle);end;
%angle= fzero(@(x) delta(x,subarea,L), 1);

% Compute the length of the new wall
length = norm(Q-P)*angle/(sin(angle));

end