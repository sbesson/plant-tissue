%% triangle
% triangle class constuctor
%%

%% Syntax   
% T = triangle(a1,a2)
%
%% Description
% Create a new object of the type triangle. This class inherits from the
% cellNetwork class
%
%% Inputs
% * a1 - the smallest angle of the triangle
% * a2 - the largest angle of the triangle
%
%% Outputs
% * T - a triangle object 
%
%% Examples
% >> T =  triangle(pi/4,pi/4);
% >> T = triangle(pi/4,pi/4,[1 1],[2 2]);
%
%
%% See also 
% * 
%
%% Author 
% Sebastien Besson.
% email address : sbesson@oeb.harvard.edu
% October 2007; Last revision: April 10, 2007

function T = triangle(varargin)

switch nargin
    case 0
        % If no input arguments, create a default object
        T.a1 = 0;
        T.a2 = 0;
        T.V1 = [0 0];
        T.V2 = [0 0];
        T = class(T,'triangle');
    case 1
        % If single argument of class triangle, return it 
        if (isa(varargin{1},'triangle')) 
            T = varargin{1};
        else
            error('Wrong argument type')  
        end 
    case 2
        % Create object using specified values
        T.a1 = varargin{1};
        T.a2 = varargin{2};
        T.V1 = [0 0];
        T.V2 = [1 0];
        T = class(T,'triangle');
    case 3
        % Create object using specified values
        a=varargin{2}-varargin{1};
        b=varargin{3}-varargin{1};
        c=varargin{3}-varargin{2};
        T.a1 = npi2pi(atan2(a(1)*b(2)-a(2)*b(1),a(1)*b(1)+a(2)*b(2)));
        T.a2 = npi2pi(atan2(c(1)*-a(2)+c(2)*a(1),-c(1)*a(1)-c(2)*a(2)));
        T.V1 = varargin{1};
        T.V2 = varargin{2};
        T = class(T,'triangle');
    case 4
        T.a1 = varargin{1};
        T.a2 = varargin{2};
        T.V1 = varargin{3};
        T.V2 = varargin{4};
        T = class(T,'triangle');
    otherwise
         error('Wrong number of arguments') 
end

end