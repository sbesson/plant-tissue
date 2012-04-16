%% findP
% Calculate the coordinates of a series of points on a curved edge.
%%

%% Syntax   
% P =findP(V,A,s);
%
%% Description 
% Calculate the coordinates of a point on a curved edge.
%
%% Inputs
% * V - a 2-by-2 matrix with the starting and ending points.
% * A - a variable containing the curvature angle of the edge with the
% cordlength.
% * s - a row vector  with the relative position of the point or points on
% the edge.  Entries of s are between 0 and 1.
%
%% Outputs
% * P - a vector with the coordinates of the points on the edge.
%
%% Examples
%  >> P =  findP([0 0;1 1], 0, 0.5);
%
%  >> P =  findP([0 0;0 1], pi/2, 0:0.1:1);
%
%% TODO
%
%% See also
% * findQ
%
%% Author  
% Jacques Dumais
% email address : jdumais@oeb.harvard.edu
% October 2007; Last revision: 31 October 2007

function [P] = findP(V,A,s)

% Check number of inputs
error(nargchk(3, 3, nargin));

% Compute the number of intermediate points and pre-allocate the P-vector
ns = length(s);

% Length of the cordlength of the edge
L = norm(V(2,:)-V(1,:));

% Unit tangent and normal vectors
T = (V(2,:) - V(1,:))/norm(V(2,:) - V(1,:)); 
N = [T(2) -T(1)];
    
if A^2 <= 0.00001 % arc approximation is valid AND required
    beta = sin(A) * L *s.*(1 - s); 
    P = ones(ns,1)*V(1,:) + L*s'*T + beta'*N;
else
    alpha = L/(2*tan(A));
    C = V(1,:) + L*T/2 - alpha*N; 
    CV = V(1,:) - C;
    P = ones(ns,1)*C + [cos(2*s'*A)*CV(1) - sin(2*s'*A)*CV(2) sin(2*s'*A)*CV(1) + cos(2*s'*A)*CV(2)] ;
end
end