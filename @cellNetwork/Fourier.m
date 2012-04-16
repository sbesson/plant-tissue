%% Fourierdescriptor
% Method of the cellNetwork class to compute the Fourier descriptor between two cellNetworks
%%

%% Syntax   
% F = Fourier(N,n)
%
%% Description
% Compute the Foturier descriptor of cells.
%
%% Inputs
% * N - a cellNetwork object
% * n - a vector of cell indices
%
%% Outputs
% * F - the result of the Fourier function
%
%% Examples
% >> F = Fourier(N,n);
%
%% See also 
% * 
%
%% Author 
% Sebastien Besson.
% email address : sbesson@oeb.harvard.edu
% July 2010; Last revision: November 21, 2008


function F = Fourier(N,n)

% Check the number of inputs
error(nargchk(1, 2, nargin));

% If no optional argument, loops over the cells
if (nargin == 1), n = 1:size(N.c,1); end;

for c = 1:length(n)
    nc = n(c);
    
    %Extracts the cell data
    C = N.c{nc};
    % Extracts the vertices and angles
    V = zeros(length(C),2);
    A = zeros(length(C),1);
    
    V(C>0,:) = N.v(N.e(C(C>0),1),:);
    A(C>0) = N.e(C(C>0),3); 
    
    V(C<0,:) = N.v(N.e(-C(C<0),2),:);
    A(C<0) = - N.e(-C(C<0),3); 
    
    dV=circshift(V,[-1 0])-V;
    L = sqrt(sum(dV.^2,2));
    smallA = A.^2<=0.00001;
    largeA = A.^2>0.00001;
    aL(smallA) = L(smallA).*(1+1/6*A(smallA).^2);
    aL(largeA) = L(largeA).*A(largeA)./sin(A(largeA));
    C=2*sin(A)./L;
    theta=atan2(dV(:,2),dV(:,1));

    dtheta = mod(circshift(theta,[-1 0])-theta,2*pi)
    %s=0:delta:1;
    s=[0 aL(1)/sum(aL) aL(1)/sum(aL) sum(aL(1:2))/sum(aL) sum(aL(1:2))/sum(aL) sum(aL(1:3))/sum(aL) sum(aL(1:3))/sum(aL)];
    theta = A(1)+[-A(1) A(1) dtheta(1)-A(2) dtheta(1)+A(2)...
        dtheta(1)+dtheta(2)-A(3) dtheta(1)+dtheta(2)+A(3) dtheta(1)+dtheta(2)+dtheta(3)-A(1)];
   %theta=theta-s*2*pi;
    
    plot(s,theta)
    set(gca,'YTick',[-pi/2 0 pi/2 pi 3*pi/2 2*pi])
end


end