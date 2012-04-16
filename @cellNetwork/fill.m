%% fill
% Fill a cellNetwork object.
%%

%% Syntax
% fill(N,n,varargin)
%
%% Description 
%  This
%
%% Inputs
% * N - a cellNetwork object
% * n - an array of cell indices
% * varargin - an optional parameter containing the degree of refinement of
% the plotted edges
%
%% Outputs
% none
%
%% Example
%
%% See also
%
%% Author
% Sebastien Besson
% email address : sbesson@oeb.harvard.edu
% November 2007; Last revision:  May 30, 2008

function h = fill(N,n,varargin)

% Loop over the cells
for c = 1:length(n)
    C= N.c{n(c)};
    celledge = [];
    for j =1:length(C)
        celledge = [celledge; edge(N,C(j))];    
    end
    % Fill the cell with the options specified by varargin
    h(c) = fill(celledge(:,1),celledge(:,2),varargin{:},'EdgeColor','none');
    hold on;
end

end
