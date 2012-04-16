%% export
% Save a cellNetwork object as an txt file.
%%
%% Syntax
% export(N,filename,varargin)
%
%% Description 
% This function exports the cellNetwork  under a txt fromat compatible with
% Surface Evolver.
%
%% Inputs
% * N - a cellNetwork structure
% * filename - the filename under which the cellNetwork will be saved
% * varargin - an optional parameter
%
%% Outputs
%
%% Example
%  >> save(N,'tissue.fe')
%
%% See also
%
%% Author
% Sebastien Besson
% email address : sbesson@oeb.harvard.edu
% March 2009; Last revision:  March 17, 2009


function export(N,filename,varargin)

fid=fopen(filename, 'wb');

if fid <= 0
  error('failed to open file');
end

fprintf(fid, '//%s\n',filename);
fprintf(fid, 'string\nspace_dimension 2\n');

fprintf(fid,'\nparameter hooke2_power  2 // the default\n');
fprintf(fid,'define edge attribute hooke_size real\n');
fprintf(fid,'quantity slinky energy method hooke2_energy global\n');

fprintf(fid, '\nvertices\n');
% Loop over vertices
for i =1:size(N.v,1)
    fprintf(fid, '%g\t%f\t%f\n',i,N.v(i,1),N.v(i,2));
end

fprintf(fid, '\nedges\n');
% Loop over edges
for i =1:size(N.e,1)
    fprintf(fid, '%g\t%g\t%g\n',i,N.e(i,1),N.e(i,2));
end


fprintf(fid, '\nfaces\n');
% Loop over edges
for i =1:length(N.c)
    fprintf(fid, '%g\t',i)
    for j =1:length(N.c{i})-1
        fprintf(fid, '%g\t',N.c{i}(j));
    end
    fprintf(fid, '%g\n',N.c{i}(length(N.c{i})))
end
edges = N.e;
edges(:,3) = 0;
newN = set(N,'Edges',edges);
fprintf(fid,'\nbodies\n');
% Loop over edges
for i =1:length(newN.c)
    if cellMoment(newN,i)>0
        fprintf(fid, '%g\t%g\tvolume %f\n',i,i,cellMoment(newN,i));
    else
        fprintf(fid, '%g\t%g\tvolume %f\n',i,-i,-cellMoment(newN,i));
    end
end

fclose(fid);
end