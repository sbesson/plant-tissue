%% subsasgn
% cellNetwork class subsref method
%%

%% Syntax   
% N = subsref(N,indices,val)
%
%% Description
% Change the property of a cellNetwork object.
%
%% Inputs
% * N - a cellNetwork object 
% * indices 
% * val - the value to be attributed corresponding property 
%
%% Outputs
% * N - the new cellNetwork object 
%
%% Examples
% >> N.vertices = zeros(1,4)
% >> N(2).v(4,1) = 1
%
%% See also 
% * subsref
%
%% Author 
% Sebastien Besson.
% email address : sbesson@oeb.harvard.edu
% June 2008; Last revision: August 4, 2009

function N = subsasgn(N,indices,val)

if length(indices) > 1
	index = indices(1);
else
	index = indices;
end

% SUBSREF 
switch index.type
    case '.'
        %val = getfield(struct(N), index.subs);
        switch index.subs
            case {'vertices','v','V'}
                if length(indices) > 1
                    x = subsasgn(N.v,indices(2:end),val);
                    N.v = x;
                else
                    N.v =val;
                end
            case {'edges','e','E'}
                if length(indices) > 1
                    x = subsasgn(N.e,indices(2:end),val);
                    N.e = x;
                else
                    N.e =val;
                end
            case {'cells','c','C'}
                if length(indices) > 1
                    x = subsasgn(N.c,indices(2:end),val);
                    N.c = x;
                else
                    N.c =val;
                end
        end

    case '()' 
        if length(indices) > 1
            x = subsasgn(N(index.subs{:}),indices(2:end),val);
            N(index.subs{:}) = x;        
        else
            N(index.subs{:}) =val;
        end
    case '{}' 
        if length(indices) > 1
            x = subsasgn(N{index.subs{:}},indices(2:end),val);
            N{index.subs{:}} = x;        
        else
            N{index.subs{:}} =val;
        end
    otherwise
        error('You do not have access to this field')
end

end
