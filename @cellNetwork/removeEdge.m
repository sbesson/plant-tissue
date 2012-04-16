%% removeEdge
% Method of the cellNetwork class to remove an edge
%%

%% Syntax   
% N = removeEdge(N,ne)
%
%% Description
% Remove an edge from a cellNetwork object.
%
%% Inputs
% * N - a cellNetwork object
% * ne - a row vector containing the indices of the edges to be removed
%
%% Outputs
% * N - a new cellNetwork object
%
%% Examples
% >> N = removeEdge(N,2); 
% remove the second edge of the network
%
%% See also 
% * addEdge
%
%% Author 
% Sebastien Besson.
% email address : sbesson@oeb.harvard.edu
% October 2007; Last revision: December 3, 2007

function N = removeEdge(N,ne)

if (length(ne) > 1)
    % Sort the list of vertices to remove in descending order
    ne = sort(ne, 'descend');
    
    % Remove vertices
    for i = 1:length(ne)
        N = removeEdge(N,ne(i));
    end

else
    % Test if the vertex exists
    if ((ne <= 0) || (ne > size(N.e,1)))
        error('removeEdge:arraySize','Edge %d does not exist',ne);
    end

    % List connected cells
    N = removeCell(N,neighbors(N,ne));
    
    % Remove the vertex
    N.e(ne,:)=[];
    
    % Shifts the edges in the remaining cells
    for i = 1:length(N.c)
        for j = 1:length(N.c{i})
            if (N.c{i}(j) > ne), N.c{i}(j) = N.c{i}(j) -1; end;
            if (N.c{i}(j) < -ne), N.c{i}(j) = N.c{i}(j) +1; end;
        end
    end

    % For debugging
    warning('removeEdge:operation','Edge %d removed',ne);

end

end