%% divide
% Divide a network according to a given solution.
%%

%% Syntax
%  N = divide(N,planes);
%
%% Description
% Divide a cell defined by its vertices and edge curvatures according to a
% Minimum vector.
%
%% Inputs
% * N - a network structure
% * planes - a solution from the divisionplanes function
%
%% Outputs
% * N - the new network structure
%% See also
% * energyLandscape
% * minima
%
%% TODO
% to be verified
%
%% Author
% Sebastien Besson
% email address: sbesson@oeb.harvard.edu
% October 2007; Last revision: October 20, 2008

function N = divide(N,planes)

% Check thenumber of inputs
error(nargchk(2, 2, nargin));

vp = zeros(1,length(planes));
ep = zeros(1,length(planes));
vq = zeros(1,length(planes));
eq = zeros(1,length(planes));

% Add the new vertices
for n = 1:length(planes)
    % Test if new vertices already exist
    [isPdefined,index] = intersect(N.v,planes(n).P,'rows');
    if isempty(isPdefined)
        N = addVertex(N,planes(n).P);
        vp(n) = size(N.v,1);
    else
        vp(n) = index;
    end
    ep(n) = N.c{planes(n).cell}(planes(n).i);
    
    % Test if new vertices already exist
    [isQdefined,index] = intersect(N.v,planes(n).Q,'rows');
    if isempty(isQdefined)
        N = addVertex(N,planes(n).Q);
        vq(n) = size(N.v,1);
    else
        vq(n) = index;
    end
    eq(n) = N.c{planes(n).cell}(planes(n).j);
end


V= [vp vq];
s =[planes(:).s1 planes(:).s2];
E = [ep eq];

% Find edges to replace
e =unique(abs(E));
% Loop over each diving edge
for n=1:length(e)
    % Find indexes that will be replaced
    i1 = find(E == e(n));
    i2 = find(E == -e(n));
    % List and sort the unique abcissas of the new vertices different from
    % 0 or 1
    [s_sort i_index] = setdiff([s(i1) 1-s(i2)],[0 1]);
    if isempty(s_sort)
        eold{n} = e(n);
        se{n} = [0 1];
    else
        sn = unique([0 s_sort 1]);
        % Compute the differences between adjacent abcissas
        ds = sn(2:end)-sn(1:end-1);
        % List the indexes to the vertices
        i =[i1 i2];
        % Define ordered list of vertices on the edge
        newV = [N.e(e(n),1) V(i(i_index)) N.e(e(n),2)];
        % Retrieve the old number of division
        division_index = N.e(e(n),4)*ones(length(ds),1);
        % Add all new edges
        N = addEdge(N,[newV(1:end-1)' newV(2:end)'  ...
            ds'*N.e(e(n),3) division_index]);
        % Save indexes of new edges
        eold{n} = size(N.e,1)-length(ds)+1:size(N.e,1);
        % Save abcissas 
        se{n} = sn;
    end
end

% Increment the order of division
division_index = (max(N.e(:,4))+1)*ones(length(vp),1);
% Add edges formed between new vertices (new division walls)
N = addEdge(N,[vq' vp' [planes(:).angle]' division_index]);
% Saves indexes to the new edges
enew =  size(N.e,1)-length(vp)+1:size(N.e,1);

% Defines daughter cells
% Daughter cell i+1:j-1 : 
for n = 1:length(planes)
    % Retrieve old data
    C = N.c{planes(n).cell};
    
    % Daughter cell 1
    % If cell is self dividing
    if planes(n).i == planes(n).j
        % Find index corresponding to the dividing edge
        ie = find(e == abs(C(planes(n).i)));
        if(C(planes(n).i)>0),
           % Find abcissas of new vertices
           [junk,index] = intersect(se{ie},[planes(n).s1 planes(n).s2]);
           % Find index of series of edges
           nij =  eold{ie}(index(1):index(2)-1);
        else
            % Find abcissas of new vertices
            [junk,index] = intersect(se{ie},[1-planes(n).s1 1-planes(n).s2]);
            % Find index of series of edges
            nij = -eold{ie}(index(2)-1:-1:index(1));
        end
        % Define daughter cell 1
        %QP->PQ
        N = addCell(N,[enew(n) nij]);
    else
        % Dividing cell excepting for new wall and edges i and j of old
        % cell
        newC1= C(planes(n).i+1:planes(n).j-1);
        % Find index corresponding to the dividing edge j
        iej = find(e == abs(C(planes(n).j)));
        if(C(planes(n).j)>0)
            % Find abcissas of new vertex
            [junk,index] = intersect(se{iej},planes(n).s2);
            % Find index of series of edges
            nj = eold{iej}(1:index-1);
        else
            [junk,index] = intersect(se{iej},1-planes(n).s2);
            nj = -eold{iej}(end:-1:index);
        end
        % Find index corresponding to the dividing edge i
        iei = find(e == abs(C(planes(n).i)));
        if(C(planes(n).i)>0) 
            % Find abcissas of new vertex
            [junk,index] = intersect(se{iei},planes(n).s1);
            % Find index of series of edges
            ni = eold{iei}(index:end);
        else
            [junk,index] = intersect(se{iei},1-planes(n).s1);
            ni = -eold{iei}(index-1:-1:1);
        end

        %i+1->j-1->s2j->QP->(1-s1)i     
        N = addCell(N,[newC1 nj enew(n) ni]);  
    end

    % Daughter cell 2
    %j+1->1->i-1->s1i->PQ->(1-s2)j
    %N = addCell(N,[C(i+1:j-1) size(N.e,1)-1 size(N.e,1)-4 size(N.e,1)-2]);
    newC2= [C(planes(n).j+1:length(C))  C(1:planes(n).i-1)];
    iej = find(e == abs(C(planes(n).j))); 
    if(C(planes(n).j)>0), 
        [junk,index] = intersect(se{iej},planes(n).s2);
        nj = eold{iej}(index:end); 
    else
        [junk,index] = intersect(se{iej},1-planes(n).s2);
        nj = -eold{iej}(index-1:-1:1);
    end

    iei = find(e == abs(C(planes(n).i)));
    if(C(planes(n).i)>0), 
        [junk,index] = intersect(se{iei},planes(n).s1);
        ni = eold{iei}(1:index-1); 
    else
        [junk,index] = intersect(se{iei},1-planes(n).s1);
        ni = -eold{iei}(end:-1:index);
    end
    N = addCell(N,[newC2 ni -enew(n) nj]);
end

% Remove mother cells
N = removeCell(N,[planes(:).cell]);

% Replace divided edges in cells
for i = 1:length(N.c)
    [ie,indexC,indexe] = intersect([N.c{i}],e);
    if ~isempty(ie)
        % Converts cell into a cell array
        temp = num2cell([N.c{i}]);
        for j =1:length(ie)
            % Replace divided edges
            temp{indexC(j)} = eold{indexe(j)};
        end
        N.c{i} = [temp{:}];
    end
    [ie,indexC,indexe] = intersect([N.c{i}],-e);
    if ~isempty(ie)
        % Converts cell into a cell array
        temp = num2cell([N.c{i}]);
        for j =1:length(ie)
            % Replace divided edges
            temp{indexC(j)} = -eold{indexe(j)}(end:-1:1);
        end
        N.c{i} = [temp{:}];
    end
end

% Clean unused vertices and edges
N = clean(N);
end