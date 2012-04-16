%% moment
% Conpute the moments of area of a given cell.
%%

%% Syntax   
% [S, C, I] = moment(V,A);
%
%% Description 
% Calculate the moments of area of a cell.
%
%% Inputs
% * V - a 2-by-2 matrix with the starting and ending points.
% * A - a variable containing the curvature angle of the edge with the
% cordlength.
%
%% Outputs
% * S - the area.
% * C - the first moment of area.
% * I - the second moment of area.
%
%% Examples
%  >> S = moment([0 0;1 1;1 0],[0 0 0])
%  >> [S, C, I] = moment([0 0;1 1;1 0],[0 0 0])
%
%% TODO
%
%% See also
% * 
%
%% Author  
% Sebastien Besson
% email address : sbesson@oeb.harvard.edu
% December 2008; Last revision: May 15, 2009

function [S, C, I] = moment(V,A)

M = V+circshift(V,[-1 0]);
dV = circshift(V,[-1 0])-V;
L = sqrt(dV(:,1).^2+dV(:,2).^2);

T = dV./repmat(L,[1 2]); % Unit  tangential vector
N = circshift(T,[0 1]); % Unit normal vector
N(:,2) = -N(:,2);

xi = V(:,1);
yi = V(:,2);
ai = xi.*circshift(yi,-1) - circshift(xi,-1).* yi;

smallA = A.^2<=0.0001;
largeA = A.^2>0.0001;
R =zeros(length(A),1);
R(largeA) = L(largeA)./abs(2*sin(A(largeA)));

%%%%%%%%%% ZERO MOMENT - AREA %%%%%%%%%%%

% Compute area of polygon
S_pol  = 1/2*sum(ai);

% Compute area of circular segments
S_cs = zeros(length(A),1);
S_cs(smallA)  = L(smallA).^2 .* A(smallA)/6;
S_cs(largeA)  = R(largeA).^2.*(A(largeA)-.5*sin(2*A(largeA)));

% Compute area of cell
S = S_pol + sum(S_cs);

if nargout == 1, return; end

%%%%%%%%%% FIRST MOMENT OF AREA %%%%%%%%%%%

% Compute first moment of area of polygon
if S_pol == 0
    C_pol  = 1/size(V,1)*sum(V,1);
else
    C_pol  = 1/(6*S_pol)*sum(M.*repmat(ai,1,2),1);
end

% Compute first moments of area of circular segments
d = zeros(length(A),1);
d(smallA) = L(smallA).*A(smallA)/10 ;
d(largeA) = R(largeA).*(4*(sin(A(largeA)).^3./(6*A(largeA)-3*sin(2*A(largeA))))...
    -cos(A(largeA)));

C_cs = M/2 + N.*repmat(d,1,2);

% Compute first moment of area of cell
C = (C_pol*S_pol+sum(C_cs.*repmat(S_cs,[1 2])))/S;

if nargout == 2, return; end

%%%%%%%%%% SECOND MOMENT OF AREA %%%%%%%%%%%

% Compute second moment of area of the polygon
IX_pol = 1/12*sum((yi.^2 + yi.*circshift(yi,-1)+ circshift(yi,-1).^2).*ai);
IY_pol = 1/12*sum((xi.^2 + xi.*circshift(xi,-1)+ circshift(xi,-1).^2).*ai);
IXY_pol =  1/24*sum((xi.*circshift(yi,-1) + 2*xi.*yi +...
    2*circshift(xi,-1).*circshift(yi,-1)+ circshift(xi,-1).*yi).*ai);
Ix_pol = IX_pol- C_pol(2)^2*S_pol;
Iy_pol = IY_pol-C_pol(1)^2*S_pol;
Ixy_pol = IXY_pol - C_pol(2)*C_pol(1)*S_pol;

% Compute second moment of area of the circular segment
D = zeros(length(A),1);
D(largeA) = 4*R(largeA).*(sin(A(largeA)).^3./(6*A(largeA)-3*sin(2*A(largeA))));
IX_cs(largeA) = R(largeA).^4/16.*(4*A(largeA)-sin(4*A(largeA))) - D(largeA).^2.*S_cs(largeA);
IY_cs(largeA) = R(largeA).^4/48.*(12*A(largeA)-8*sin(2*A(largeA))+sin(4*A(largeA)));
IX_cs(smallA) = L(smallA).^4.*A(smallA).^3/1400;
IY_cs(smallA) = L(smallA).^4.*A(smallA)/120;

phi = atan2(-dV(:,2),-dV(:,1))';

Ix_cs = (IX_cs+IY_cs)/2+(IX_cs-IY_cs)/2.*cos(2*phi);
Iy_cs = (IX_cs+IY_cs)/2-(IX_cs-IY_cs)/2.*cos(2*phi);
Ixy_cs = (IX_cs-IY_cs)/2.*sin(2*phi);

% Parallel axis theorem for the polygon
Tx_pol = (C(2)-C_pol(2))^2*S_pol;
Ty_pol = (C(1)-C_pol(1))^2*S_pol;
Txy_pol =(C(2)-C_pol(2))*(C(1)-C_pol(1))*S_pol;

% Parallel axis theorem for the circular segment
Tx_cs = (C(2)-C_cs(:,2)).^2.*S_cs;
Ty_cs = (C(1)-C_cs(:,1)).^2.*S_cs;
Txy_cs =(C(2)-C_cs(:,2)).*(C(1)-C_cs(:,1)).*S_cs;

Ix = Ix_pol + Tx_pol + sum(Ix_cs) + sum(Tx_cs);% + sum(T_cs,3);
Iy = Iy_pol + Ty_pol + sum(Iy_cs) + sum(Ty_cs);% + sum(T_cs,3);

%Iy_pol +  sum(Iy_cs)
Ixy = Ixy_pol + Txy_pol + sum(Ixy_cs) + sum(Txy_cs);% + sum(T_cs,3);
I = [Ix -Ixy; -Ixy Iy];
end