%% removeVertex
% Method of the cellNetwork class to remove a vertex
%%

%% Syntax   
% N = removeVertex(N,nv)
%
%% Description
% Remove a vertex from a cellNetwork object.
%
%% Inputs
% * N - a cellNetwork object
% * nv - a row vector containing the indices of the vertices to be removed
%
%% Outputs
% * N - a new cellNetwork object
%
%% Examples
% >> N = removeVertex(N,1); 
% remove the first vertex of the network
%
%% See also 
% * addVertex
%
%% Author 
% Sebastien Besson.
% email address : sbesson@oeb.harvard.edu
% October 2007; Last revision: January 25, 2007

function N = removeVertex(N,nv)


if (length(nv) > 1)
    % Sort the list of vertices to remove in descending order
    nv = sort(nv, 'descend');
    
    % Remove vertices
    for i = 1:length(nv)
        N = removeVertex(N,nv(i));
    end

else
    % Test if the vertex exists
    if ((nv <= 0) || (nv > size(N.v,1)))
        error('removeVertex:arraySize','Vertex %d does not exist',nv);
    end

    % List connected edges
    e_index = find(transpose(N.e(:,1:2) == nv));
    e_index = round(e_index/2);

    % Remove connected edges
    if (~isempty(e_index)), N = removeEdge(N,e_index); end;
    
    % Remove the vertex
    N.v(nv,:)=[];

    % Shifts the vertices in theremaining edges
    for i = 1:size(N.e,1)
        if (N.e(i,1) > nv), N.e(i,1) = N.e(i,1) -1; end;
        if (N.e(i,2) > nv), N.e(i,2) = N.e(i,2) -1; end;
    end

    % For debugging
    warning('removeVertex:operation','Vertex %d removed',nv);

end

end