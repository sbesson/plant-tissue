%% subsref
% cellNetwork class subsref method
%%

%% Syntax   
% val = subsref(N,indices)
%
%% Description
% Return the property of a cellNetowrk object.
%
%% Inputs
% * N - a cellNetwork object 
% * indices 
%
%% Outputs
% * val - the value of the corresponding property 
%
%% Examples
% >> N.vertices
% >> N(2).v(4)
%
%% See also 
% * subsasgn
%
%% Author 
% Sebastien Besson.
% email address : sbesson@oeb.harvard.edu
% June 2008; Last revision: August 4, 2009

function val = subsref(N,indices)

if length(indices) > 1
	index = indices(1);
else
	index = indices;
end
% SUBSREF 
switch index.type
    case '.'
        switch index.subs
            case {'vertices','v','V'}
                val = N.v;
            case {'edges','e','E'}
                val = N.e;
            case {'cells','c','C'}
                val = N.c;
        end
    case '()' 
        val = N(index.subs{:});
    case '{}' 
        val = N{index.subs{:}};
    
    otherwise
        error('You do not have access to this field')
end

if length(indices) > 1
	val = subsref(val, indices(2:end));
end
end
