%% col2cell
% Converts a coleochaete object into a cellNetwork object.
%%

%% Syntax   
% N = col2cell(C)
%
%% Description
% Create a new cellNetwork out of a coleoachete object
%
%% Inputs
% * C - a coleochaete object
%
%% Outputs
% * N - a cellNetwork object 
%
%% Examples
% >> N = col2cell(coleochaete(.2,pi/4));
% >> N = col2cell(coleochaete(.2,pi/4,[0 0],[1 0]));
%
%
%% See also 
% * cell2col 
%
%% Author 
% Sebastien Besson.
% email address : sbesson@oeb.harvard.edu
% January 2008; Last revision: July 7, 2009

function N = col2cell(C)

theta= C.theta;
h= C.h;
D = norm(C.V1-C.V2);

if D == 0 
    R = sqrt(2/(theta*h*(2-h))); % not precised - scaled to unit area
    N = [0 1];
    T= [1 0];
    d = 4/3*R*sin(theta/2)/theta*(1+(1-h)^2/(2-h));
    P1 =  [R*cos(theta/2)-d -R*sin(theta/2)];
    P2 =  [R*cos(theta/2)-d R*sin(theta/2)];
else
    R = norm(C.V1-C.V2)/(2*sin(theta/2));
    N = (C.V2-C.V1)/norm(C.V1-C.V2);
    T = [N(2) -N(1)];
    P1 = C.V1;
    P2 = C.V2;
end
    
if h==1
    V = [P1;P2;P1-h*R*cos(theta/2)*T+h*R*sin(theta/2)*N];
    E = [1 2 theta/2 1;
        2 3 0 3;
        3 1 0 3];
    C ={[1 2 3]};
else
    V = [P1;P2;
            P2-h*R*cos(theta/2)*T-h*R*sin(theta/2)*N;
            P1-h*R*cos(theta/2)*T+h*R*sin(theta/2)*N;];
    E = [1 2 theta/2 1;
        2 3 0 3;
        3 4 -theta/2 2;
        4 1 0 3];
    C ={[1 2 3 4]};
end
N =  cellNetwork(V,E,C);
end