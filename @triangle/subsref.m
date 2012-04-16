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

function val = subsref(T,indices)

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
            case 'a1'
                val = T.a1;
            case 'a2'
                val = T.a2;
            case 'V1'
                val = T.V1;
            case 'V2'
                val = T.V2;
        end
        %if length(index) > 1 &&  strcmp(index(2).type,'()')
        %    val = val(index(2).subs{:});
        %end 

    case '()' 
        %if (length([index.subs]) ~= 1)
        %    error('Only single indexing is supported.');
        %end      
        val = T(index.subs{:});
    case '{}' 
        %if (length([index.subs]) ~= 1)
        %    error('Only single indexing is supported.');
        %end      
        val = T{index.subs{:}};
    
    otherwise
        error('You do not have access to this field')
     
end

if length(indices) > 1
	val = subsref(val, indices(2:end));
end

end
