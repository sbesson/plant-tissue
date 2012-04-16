%% strands
% Determines the strands connecting the nucleus to the mother wall of a cell.
%%

%% Syntax
% S = strands(N,n)
% 
%% Description
% Compute the
%
%% Inputs
% * N - a cellNetwork object
% * n - an array of cell indices
%
%% Outputs
% * 
%
%% Example
%	>> strands = strands(N,1); 
%
%% See also
% * cellCentroid
%
%% Author
% Sébastien Besson
% email address : jdumais@oeb.harvard.edu
% May 2008; Last revision: May 26, 2008

function [S cC] = strands(N,n)

% Check the number of inputs
error(nargchk(1, 2, nargin));

% If no optional argument, loop over all cells
if (nargin == 1), n = 1:length(N.c); end;

cC = zeros(length(n),2);
%planes = struct('cell',0,'i',[],'s1',[],'P',[],'j',[],'s2',[],'Q',[],'walllength',[]);
%strands = struct('cell',0,'i',[],'s1',[],'P',[],'j',[],'s2',[],'Q',[],'walllength',[]);

% Loop over the different cells
for c = 1:length(n)
    % Extract the cell data
    C = N.c{n(c)};
    
    % Extract the vertices and angles
    V = zeros(length(C),2);
    A = zeros(length(C),1);
    
    index1 = logical(C>0);
    V(index1,:) = N.v(N.e(C(index1),1),:);
    A(index1) = N.e(C(index1),3); 
    
    index2 = logical(C<0);
    V(index2,:) = N.v(N.e(-C(index2),2),:);
    A(index2) = - N.e(-C(index2),3); 

    % Add a supplementary vertex to easy the writing of the program
    T = circshift(V,[-1 0])-V;
    L = sqrt(T(:,1).^2+T(:,2).^2);
    N= circshift(T,[0 1]);
    N(:,2) = -N(:,2);

    aL = zeros(size(V,1),1);
    h = zeros(size(V,1),1);
    % Compute area of circular segments
    ind1 = A.^2<=0.00001;
    aL(ind1) = L(ind1).*(1+1/6*A(ind1).^2);
    h(ind1) = A(ind1)/6;

    ind2 = A.^2>0.00001;
    aL(ind2) = L(ind2).*A(ind2)./sin(A(ind2));
    h(ind2) = 1./A(ind2) -1./tan(A(ind2));

    points = V+T/2+ repmat(h,1,2).*N/2;
    cC(c,:) = sum(repmat(aL,1,2).*points,1)/sum(aL);

    
    cCV = V-repmat(cC(c,:),size(V,1),1);
    s = zeros(size(V,1),1);
    s(ind1) = 1/2+1/2*(2*dot(cCV(ind1,:),T(ind1,:),2)+L(ind1).^2)./...
         (2*A(ind1).*dot(cCV(ind1,:),N(ind1,:),2)- (L(ind1).^2));

    s(ind2) = 1/2+1./(2*A(ind2)).*atan((2*dot(cCV(ind2,:),T(ind2,:),2)+L(ind2).^2)./...
        (2*dot(cCV(ind2,:),N(ind2,:),2)- (L(ind2).^2)./tan(A(ind2))));
    
    % If found point is outside the arc of circle
    outofbounds = (s<0) | (s>1);
    D = sqrt(cCV(:,1).^2+cCV(:,2).^2);
    delta = D - circshift(D,[-1 0]); % CV1-CV2
    % if delta>0, s=1 else s=0
    s(outofbounds) = delta(outofbounds)>0; 
    
    P=zeros(size(V));
    V2 = circshift(V,[-1 0]);
    for i = 1:size(V,1)
        P(i,:) = findP([V(i,:);V2(i,:)],A(i),s(i));
    end
    
    S(c).i=1:size(V,1);
    S(c).P = P;
    S(c).s = s;
    S(c).centroid = cC;
    S(c).weight = aL/sum(aL);
end