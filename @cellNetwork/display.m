%% display
% Method of the cellNetwork class to display information
%%

%% Syntax   
% l = display(N)
%
%% Description
% Print information about the subelements of the network.
%
%% Inputs
% * N - a cellNetwork object
%
%% Outputs
% * l - a string
%
%% Examples
% >> l = display(N); 
% print information about the N cellNetwork object.
%
%% See also 
% * addEdge
%
%% Author 
% Sebastien Besson.
% email address : sbesson@oeb.harvard.edu
% October 2007; Last revision: August 12, 2008

function l = display(N)

N=N(1); % Display only the first element of an array.

l=sprintf('Vertices\n');
if (size(N.v,1))  
    for i = 1:size(N.v,1)
        l = sprintf('%sVertex %d: %f %f\n',l,i,N.v(i,1),N.v(i,2));
    end
else
    l = [l 'No vertex'];
end
    
l=[l sprintf('\nEdges\n')];
if (size(N.e,1))
    for i = 1:size(N.e,1)
        l = sprintf('%sEdge %d: from vertex %d to %d with an angle %f\n',...
            l,i,N.e(i,1),N.e(i,2),N.e(i,3));
    end
else
      l = [l 'No edge'];    
end

l=[l sprintf('\nCells\n')];

if ~isempty(N.c)
    for i = 1:length(N.c)
        l = sprintf('%sCell %d: ',l,i);
        for j = 1:length(N.c{i})
            if (N.c{i}(j)>0)
                l = sprintf('%s%d ->',l,N.e(N.c{i}(j),1));
            else
                l = sprintf('%s%d ->',l,N.e(-N.c{i}(j),2));
            end
        end
        
        if (N.c{i}(j)>0)
            l = sprintf('%s%d\n',l,N.e(N.c{i}(j),2));
        else
            l = sprintf('%s%d\n',l,N.e(-N.c{i}(j),1));
        end

    end
else
        l = [l 'No cell'];
end

disp(l);
end