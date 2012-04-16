%% eq
% Method of the cellNetwork class to compare two cellNetworks
%%

%% Syntax   
% test = eq(N1,N2)
%
%% Description
% Add a new cell to a cellNetwork object.
%
%% Inputs
% * N1 - a cellNetwork object
% * N2 - a cellNetwork object
%
%% Outputs
% * test - the result of the equal function
%
%% Examples
% >> eq(N1,N2);
% >> N1 == N2;
%
%% See also 
% * 
%
%% Author 
% Sebastien Besson.
% email address : sbesson@oeb.harvard.edu
% November 2007; Last revision: November 21, 2008

function test = eq(N1,N2)

% Check the number of inputs
error(nargchk(2, 2, nargin));

if ~isa(N1,'cellNetwork') && ~isa(N2,'cellNetwork')
    error('cellNetwork:distance',...
        'Both arguments must be cellNetwork objects.') 
end

if ~all(size(N1.v) == size(N2.v)), test = 0; return; end
if ~all(size(N1.e) == size(N2.e)), test = 0; return; end
C1=[N1.c{:}];
C2=[N2.c{:}];
if length(C1) ~= length(C2), test = 0; return; end

test = (all(N1.v==N2.v) & all(N1.e==N2.e) & all(C1==C2));

end