%% distance
% Method of the cellNetwork class to compute the distance between two cellNetworks
%%

%% Syntax   
% D = distance(N1,N2)
%
%% Description
% Compute the Hausdorff distance between two cellNetworks.
%
%% Inputs
% * N1 - a cellNetwork object
% * N2 - a cellNetwork object
%
%% Outputs
% * D - the result of the distance function
%
%% Examples
% >> D = distance(N1,N2);
%
%% See also 
% * 
%
%% Author 
% Sebastien Besson.
% email address : sbesson@oeb.harvard.edu
% November 2007; Last revision: November 21, 2008

%%%
%%% HAUSDORFF: Compute the Hausdorff distance between two point clusters 
%%%            in an arbitrary dimensional vector space.
%%%            H(A,B) = max(h(A,B),h(B,A)), where
%%%            h(A,B) = max(min(d(a,b))), for all a in A, b in B,
%%%            where d(a,b) is a L2 norm.
%%%   dist = hausdorff( A, B )

function dist = distance(N1,N2)

% Check the number of inputs
error(nargchk(2, 2, nargin));

if ~isa(N1,'cellNetwork') && ~isa(N2,'cellNetwork')
    error('cellNetwork:distance',...
        'Both arguments must be cellNetwork objects.') 
end

dist = max(compute_dist(N1,N2),compute_dist(N2,N1));
%dist = dist/abs(cellArea(N1));

    %%% Compute distance 
    function distance = compute_dist(N1,N2)
        
        ages = unique(N1.e(:,4));
        d = zeros(length(ages),1);
        for i =1:length(ages)
            
            res =arrayfun(@(x) edge(N1,x),find(N1.e(:,4)==ages(i))','UniformOutput',false);
            B1 = cell2mat(res');  

            res =arrayfun(@(x) edge(N2,x),find(N2.e(:,4)==ages(i))','UniformOutput',false);
            B2 = cell2mat(res');  
            distances = zeros(length(B2),1);

            for k = 1 : length(B1)
                B12 = ones(length(B2),1)*B1(k,:);
                D = (B12-B2).^2;
                D = sqrt(D*ones(2,1));
                distances(k) = min(D);
            end
            d(i) = max(distances);
        end
        
        distance = sum(d);
    end

end