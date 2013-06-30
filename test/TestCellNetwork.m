% Unit tests for the cellNetwork class
% Require MATLAB xUnit Test Framework to be installed
% http://www.mathworks.com/matlabcentral/fileexchange/22846-matlab-xunit-test-framework


classdef TestCellNetwork < TestCase
    
    properties
        v
        e
        c
        tissue
    end
    
    methods
        function self = TestCellNetwork(name)
            self = self@TestCase(name);
        end
        
        % Test constructors
        function checkNetwork(self)
            assertTrue(isa(self.tissue, 'cellNetwork'));
            assertEqual(self.tissue.v, self.v);
            assertEqual(self.tissue.e, self.e);
            assertEqual(self.tissue.c, self.c);
        end
        
        function createMonoCellNetwork(self, nv)
            self.v = zeros(nv, 2);
            self.e = [(1:nv)' circshift((1:nv)',-1) zeros(nv,1) ones(nv,1)];
            self.c = {1:nv};
            self.tissue = cellNetwork(self.v, self.e, self.c);
        end
        
        function testEmptyNetwork(self)
            self.tissue = cellNetwork();
            self.v = [];
            self.e = [];
            self.c = {};
            self.checkNetwork();
        end
        
        
        function testTriangleNetwork(self)
            self.createMonoCellNetwork(3);
            self.checkNetwork();
        end
        
        function testSquareNetwork(self)
            self.createMonoCellNetwork(4);
            self.checkNetwork();
        end
        
        function testSquareNetworkCopy(self)
            self.createMonoCellNetwork(4);
            self.tissue = cellNetwork(self.tissue);
            self.checkNetwork();
        end
    end
end
