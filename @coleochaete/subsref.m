%% subsred
% microscope class subsref method
%%

%% Syntax   
% val = subsref(m,field)
%
%% Description
% Return the property of a microscope object.
%
%% Inputs
% * m - a microscope object 
% * field - a string containing the field
%
%% Outputs
% * val - the corresponding property 
%
%% Examples
% >> m.dimensions
% >> roi = m.roi
%
%% See also 
% * get
%
%% Author 
% Sebastien Besson.
% email address : sbesson@oeb.harvard.edu
% June 2008; Last revision: October 20, 2008

function val = subsref(C,indices)

if length(indices) > 1
	index = indices(1);
else
	index = indices;
end
% SUBSREF 
switch index.type
    case '.'
        switch index.subs
            case 'h'
                val = C.h;
            case 'theta'
                val = C.theta;
            case 'V1'
                val = C.V1;
            case 'V2'
                val = C.V2;
        end
    case '()' 
        if (length(index.subs) ~= 1)
            error('Only single indexing is supported.');
        end
        val = C(index.subs{:});
    case '{}' 
        val = C{index.subs{:}};
    otherwise
        error('You do not have access to this field')    
end

if length(indices) > 1
	val = subsref(val, indices(2:end));
end
end