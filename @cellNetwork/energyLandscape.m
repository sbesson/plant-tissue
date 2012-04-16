%% energyLandscape
% Compute the energy landscape of a dividing cell.
%%

%% Syntax
% planes = energyLandscape(N, n,delta)
%
%% Description
% The boundary of the cell is discretized with a given resolution. For each
% point on the boundary, dividing solutions are probed on each wall. If
% such a solution exists, it is stored in a structure.
% 
%% Inputs
% * N - a cellNetwork object
% * n - an array of indices
% * delta - an optional input containing the resolution for the
% discretisation of the cell.
%
%% Outputs
% planes is an array of structures containing the following information:
% 
% * arcposition is the position of the point on the length of the outer
% wall.
% * i is the index of the first wall.
% * j is the index of the second wall.
% * P are the coordinates of the extrmity on the first wall.
% * Q are the coordinates of the extremity on the second wall. 
% * angle is the angle formed between the wall and the cord.
% * length is the length of the new wall.
%
%% Example
% >> Solutions =  energyLandscape([0 0;0 1; 1 1; 1 0], [0 0 0 0]);
%
% >> Solutions =  energyLandscape([0 0;0 1; 1 1; 1 0], [0 0 0 0],.01);
%
%% See also
% * energyAnimation
% * plotLandscape
%
%% Author
% Jacques Dumais, Sebastien Besson
% email address : sbesson@oeb.harvard.edu
% October 2007; Last revision: 31 October 2007

function planes = energyLandscape(N,n,delta)

% Check the number of arguments
error(nargchk(2, 3, nargin))

% If no optional argument, loops over the cells
if (nargin == 1), n = 1:size(N.c,1); end;

% If empty second argument
if isempty(n), planes = []; return; end;

% If delta is not precised, attribute the default value to the
% discretization resolution
if (nargin ==2) , delta = 0.005; end

%minArray = struct('cell',0,'i',0,'j',0,'s1',0.0,'s2',0.0,'P',[],'Q',[],'angle',0.0,'walllength',0.0);
planes = struct('cell',0,'i',0,'j',0,'s1',0.0,'s2',0.0,'arcposition',0.0,...
    'P',[],'Q',[],'angle',0.0,'length',0.0,'color',0.0);

for c = 1:length(n)
    nc = n(c);
    
    %Extracts the cell data
    C = N.c{nc};
    % Extracts the vertices and angles
    V = zeros(length(C),2);
    A = zeros(length(C),1);
    
    for i = 1:length(C)
        if (C(i)>0)
            V(i,:) = N.v(N.e(C(i),1),:);
            A(i) = N.e(C(i),3); 
            continue; 
        else
            V(i,:) = N.v(N.e(-C(i),2),:);
            A(i) = - N.e(-C(i),3); 
            continue; 
        end
    end

    % Compute area of the cell
    cArea = area(V,A);
    
    % Initialize the solutions counter
    count = 0;    
    
    % Add a supplementary vertex to easy the writing of the program
    V = [V;V(1,:)];

    cL = zeros(size(V,1)-1,1); % chord length between points
    aL = zeros(size(V,1)-1,1); % arc length between points
    np = zeros(size(V,1)-1,1); % number of points discretizing each edge


    for i = 1:size(V,1)-1
        cL(i) = norm(V(i+1,:)-V(i,:)); 
        if A(i)^2 <= 0.00001 % arc approximation is valid AND required
           aL(i) = norm(V(i+1,:)-V(i,:));
        else
           aL(i) = norm(V(i+1,:)-V(i,:))*A(i)/sin(A(i)); 
        end

        np(i) = round(aL(i)/delta);  % number of points to discretize the edge
    end   

    %n = sum(np); % total number of points making up the cell outline    
    %arcpositions = arcpositions/arcpositions(n); c


    % Loop on all the couple of edges (i,j)
 
    for i = 1:size(V,1)-1
        for s0 = 0:delta:1
            for j = 1:size(V,1)-1      
                if i==1 
                    arcposition = s0*aL(1)/sum(aL(:));  
                else
                    arcposition = (s0*aL(i) + sum(aL(1:i-1)))/sum(aL(:));
                end

                if (i == j)
                    s = fminsearch(@(s) findLength2(s0,s,V,A,i,j,cArea/2),(1-s0/2));
                else
                    s = fminsearch(@(s) findLength2(s0,s,V,A,i,j,cArea/2),.5);
                end

                % Validity tests

                % Test the position of the points P and Q
                %isValid1 = (s >= 0 && s <= 1);
                isValid1 = (s >= -1e-4 && s <= 1+1e-4);
                if ~isValid1, continue; end;

                [walllength,P,Q,angle] = findLength2(s0,s,V,A,i,j,cArea/2);

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
                planes(c,count).s1 = s0;
                planes(c,count).j = j;
                planes(c,count).s2 = s;
                planes(c,count).arcposition = arcposition;
                planes(c,count).P = P;
                planes(c,count).Q = Q;
                planes(c,count).angle = angle;
                planes(c,count).length = walllength;  
                %plot([planes(:).arcposition],[planes(:).length],'o');
                drawnow
            end
        end
    end
    
    % Sort the results according to the length of the wall
    %[junk,index] = sort([planes(c,1:count).walllength]);
    %minArray(c,1:count) = minArray(c,index);
end


