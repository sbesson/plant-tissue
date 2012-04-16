%% divisionplanes
% Compute the division planes for of a dividing cell.
%%
%% Syntax  
% planes = divisionplanes(N,n,area_ratio);
%
%% Description 
% For each couple of edges (i,j), the length of a wall dividing the cell
% into two daughter cells of equal area is minimized. If the minimum does
% not start or end at a vertex, it is stored in a solution array. 
%
%% Inputs
% * N - a cellNetwork object
% * n - an array of indices
% * area_ratio - an area ratio
%
%% Outputs
% planes is an array of structures containing the following information:
%
% * [i,j] are the edges.
% * [Q,P] are the values of the vertices defining the new wall.
% * angle is the value of the angle of the new angle defined for the ede
% [QP]
% * [s1,s2] are the position of the 
%
%% Examples
% >> divisionplanes([0 0; 1 0; 0 1], [0 pi/4 0]);
% returns three minimal solutions corresponding to the division of a
% quadrant cell.
%
%% See also 
% * cellArea
% * energyLandscape
% * findLength
%
%% Author
% Sebastien Besson
% email address: sbesson@oeb.harvard.edu
% October 2007; Last revision: July 22, 2010

function planes = divisionplanes(N,n,ratio)

% Check the number of inputs
error(nargchk(1, 3, nargin));

% If no optional argument, loops over the cells
if (nargin < 3), ratio = .5; end;
if (nargin < 2), n = 1:length(N.c); end;

% If empty second argument
 if isempty(n), planes = []; return; end;

planes = struct('cell',0,'i',0,'j',0,'s1',0.0,'s2',0.0,'P',[],'Q',[],...
    'angle',0.0,'walllength',0.0,'probability',0.0);

% Probability distribution
alpha = 20.6;

for c = 1:length(n)
    nc = n(c);
    
    %Extracts the cell data
    C = N.c{nc};
    % Extracts the vertices and angles
    V = zeros(length(C),2);
    A = zeros(length(C),1);
    
    index1 = logical(C>0);
    V(index1,:) = N.v(N.e(C(index1),1),:);
    A(index1) = N.e(C(index1),3); 
    
    index2 = logical(C<0);
    V(index2,:) = N.v(N.e(-C(index2),2),:);
    A(index2) = - N.e(-C(index2),3); 
    
    % Compute area of the cell
    cArea = area(V,A);

    % Initialize the solutions counter
    count = 0;

    % Add a supplementary vertex to easy the writing of the program
    V = [V;V(1,:)];

    % Loop on all the couple of edges (i,j)
    for i = 1:size(V,1)-1
        for j = i:size(V,1)-1
           if (i == j)
               s = fminsearch(@(s) findLength(s,V,A,i,j,ratio*cArea),[0.25 0.75]);
           else
               s = fminsearch(@(s) findLength(s,V,A,i,j,ratio*cArea),[0.5 0.5]);
           end

           % Validity tests
           % Test the position of the points P and Q
           % Allow s to be equal to 0
           tol =1e-3;
           isValid1 = (s(1) >= -tol && s(1) < 1 && s(2) >= -tol && s(2) <1);
           if ~isValid1, continue; end;

           if s(1)<=tol, s(1) =0;end
           if s(2)<=tol, s(2) =0;end
           [walllength,P,Q,angle] = findLength(s,V,A,i,j,ratio*cArea);

           % Solution must be an arc of a circle
           % Angle must belong to the interval ]-pi;pi[
           isValid2 = (angle < pi && angle > -pi);
           if ~isValid2, continue; end;

           % Test if the new wall does not intersect other walls.
           %isValid3 = intersection(P,Q,angle,V,A);
           %if ~isValid3, continue; end;

           % Solutions are stored in the array
           count = count +1;
           planes(c,count).cell = nc;
           planes(c,count).i = i;
           planes(c,count).s1 = s(1);
           planes(c,count).j = j;
           planes(c,count).s2 = s(2);
           planes(c,count).P = P;
           planes(c,count).Q = Q;
           planes(c,count).angle = angle;
           planes(c,count).walllength = walllength; 
        end      
    end
         
    Z = sum(exp(-alpha*[planes(c,1:count).walllength]/sqrt(cArea)));
    for j=1:count
        planes(c,j).probability = exp(-alpha*planes(c,j).walllength/sqrt(cArea))/Z;
    end
    % Sort the results according to the length of the wall
    [junk,index] = sort([planes(c,:).walllength]);
    planes(c,1:length(index)) = planes(c,index);
end

end
