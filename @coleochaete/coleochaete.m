%% coleochaete
% coleochaete class constuctor
%%

%% Syntax   
% C = coleochaete(x,theta)
%
%% Description
% Create a new object of the type coleochaete.
%
%% Inputs
% * x - the ratio of the length over the total radius
% * theta - the angle of the cell
%
%% Outputs
% * C - a coleochaete object 
%
%% Examples
% >> C =  coleochaete(.1,pi/4);
%
%% See also 
% * triangle
%
%% Author 
% Sebastien Besson.
% email address : sbesson@oeb.harvard.edu
% November 2008; Last revision: February 17, 2008

function C = coleochaete(varargin)

switch nargin
    case 0
        % If no input arguments, create a default object
        C.h = 0.5;
        C.theta = pi/4;
        C.V1 = [0 0];
        C.V2 = [0 0];
        C = class(C,'coleochaete');
    case 1
        % If single argument of class triangle, return it 
        if (isa(varargin{1},'coleochaete')) 
            C = varargin{1};
        else
            error('Wrong argument type')  
        end 
    case 2
        % Create object using specified values
        C.h = varargin{1};
        C.theta = varargin{2};
        C.V1 = [0 0];
        C.V2 = [0 0];
        C = class(C,'coleochaete');
    case 4
        C.h = varargin{1};
        C.theta = varargin{2};
        C.V1 = varargin{3};
        C.V2 = varargin{4};
        C = class(C,'coleochaete');

    otherwise
         error('Wrong number of arguments') 
end

end