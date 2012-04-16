%% centroid
% Conpute of the center of mass of a given cell.
%%

%% Syntax   
% C = centroid(V,A);
%
%% Description 
% Calculate the coordinate of the centroid of a cell.
%
%% Inputs
% * V - a 2-by-2 matrix with the starting and ending points.
% * A - a variable containing the curvature angle of the edge with the
% cordlength.
%
%% Outputs
% * C - the coordinates of the cell centroid.
%
%% Examples
%  >> C = centroid([0 0;1 1;1 0],[0 0 0])
%
%% TODO
%
%% See also
% * 
%
%% Author  
% Sebastien Besson
% email address : sbesson@oeb.harvard.edu
% October 2007; Last revision: December 22, 2008

function [Cmin P] = centroid2(V,A)
V2 = circshift(V,[-1 0]);
T = circshift(V,[-1 0])-V;
L = sqrt(T(:,1).^2+T(:,2).^2);
T = T./repmat(L,1,2);
N= circshift(T,[0 1]);
N(:,2) = -N(:,2);

aL = zeros(size(V,1),1);
d = zeros(size(V,1),1);

% Compute area of circular segments
ind1 = A.^2<=0.00001;
aL(ind1) = L(ind1).*(1+1/6*A(ind1).^2);

ind2 = A.^2>0.00001;
aL(ind2) = L(ind2).*A(ind2)./sin(A(ind2));


for i = 1:size(V,1)
    P(i,:) = findP([V(i,:);V2(i,:)],A(i),1/2);
end
Cmin = sum(repmat(aL,1,2).*P,1)/sum(aL);

%     function [E P] = energy(C)
%         CV = V-repmat(C,size(V,1),1);
%         CV_T = dot(CV,T,2);
%         CV_N = dot(CV,N,2);
%         s(ind1) = 1/2+1/2*(2*CV_T(ind1)+L(ind1))./(2*A(ind1).*CV_N(ind1)-L(ind1));
%         s(ind2) = 1/2+1./(2*A(ind2)).*atan(tan(A(ind2)).*(2*CV_T(ind2)+L(ind2))./...
%             (2*tan(A(ind2)).*CV_N(ind2)-L(ind2)));
%         s =  max(min(s,1),0);
%         for i = 1:size(V,1)
%             P(i,:) = findP([V(i,:);V2(i,:)],A(i),s(i));
%             d(i) = (C(1)-P(i,1))^2 +(C(2)-P(i,2))^2;
%         end
%         E = sum(aL.*d);
%     end

%     function [E C P] = energy(s)
%          s =  max(min(s,1),0);
%         for i = 1:size(V,1)
%             P(i,:) = findP([V(i,:);V2(i,:)],A(i),s(i));
%         end
%         C = sum(repmat(aL,1,2).*P,1)/sum(aL);
% 
%         CP = P-repmat(C,size(P,1),1);
%         d = CP(:,1).^2+CP(:,2).^2;
%         E = sum(aL.*d);
%     end
% 
% 
% smin = fminsearch(@(x) energy(x),.5*ones(size(V,1),1));
% [Emin Cmin, P] = energy(smin);
%energy([0.5 .2]);

end