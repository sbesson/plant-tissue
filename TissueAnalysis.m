function varargout = TissueAnalysis(varargin)
% TissueAnalysis M-file for TissueAnalysis.fig
%      TissueAnalysis, by itself, creates a new TissueAnalysis or raises the existing
%      singleton*.
%
%      H = TissueAnalysis returns the handle to a new TissueAnalysis or the handle to
%      the existing singleton*.
%
%      TissueAnalysis('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TissueAnalysis.M with the given input arguments.
%
%      TissueAnalysis('Property','Value',...) creates a new TissueAnalysis or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TissueAnalysis_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TissueAnalysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% Cells the above text to modify the response to help TissueAnalysis

% Last Modified by GUIDE v2.5 15-Apr-2012 16:03:43

% See also: GUIDE, GUIDATA, GUIHANDLES

% Begin initialization code - DO NOT CELLS
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TissueAnalysis_OpeningFcn, ...
                   'gui_OutputFcn',  @TissueAnalysis_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

% End initialization code - DO NOT CELLS

% --- Executes just before TissueAnalysis is made visible.
function TissueAnalysis_OpeningFcn(hObject, eventdata, handles)

clc
warning off all
handles.pathname = '';
handles.filename = 0;
handles.tissue = cellNetwork();
handles.tissueplot = plot(handles.tissue);
handles.selectedvertex = 0;
handles.selectededge = 0;
handles.selectedcell = 0;
handles.initialpoint = [0 0];
handles.firstcellvertex = 0;
handles.newcell = [];
handles.buttons = [handles.create handles.edit handles.clear handles.defineCell handles.autodefineCell];

% Default marker and line properties
handles.defaultmarkersize = 4;
handles.selectedmarkersize = 10;
handles.defaultlinewidth = 1;
handles.selectedlinewidth = 4;
set(handles.agelabel,'BackgroundColor',coloredge(1));

handles.output = hObject; % Default output

guidata(hObject, handles);  % Update handles structure

% --- Outputs from this function are returned to the command line.
function varargout = TissueAnalysis_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

% --------------------------------------------------------------------
function OpenImage_Callback(hObject, eventdata, handles)

set(handles.toggle_openimage,'Value',0);    % Set the toggle button on

[filename,pathname] = uigetfile({'*.bmp;*.tif;*.png;*.jpg'},...
    'Choose image file to treat',handles.pathname);

% If valid filename launch the openimage function
if ~isequal(filename,0), 
    handles.filename = filename;
    handles.pathname = pathname;
    guidata(hObject, handles);  % Update handles structure
    openimage(hObject,eventdata); 
end

% --------------------------------------------------------------------
function ImportData_Callback(hObject, eventdata, handles)

[filename,path] = uigetfile({'*.mat'},'Choose data file to import',handles.pathname);

% If valid filename launch the openimage function
if ~isequal(handles.filename,0)
    opendata(hObject,eventdata,[path filename]); 
end

function openimage(hObject,eventdata)

handles = guidata(hObject); % Retrieve handles structure

% Reinitialize cellNetwork and graphic objects
delete([handles.tissueplot.edges handles.tissueplot.vertices]);
handles.tissue = cellNetwork();
handles.tissueplot = plot(handles.tissue);
handles.selectedvertex = 0;
handles.selectededge = 0;
handles.selectedcell = 0;

set(handles.buttons,'Enable','on','Value',0);  % Initialize toggle buttons

% Read the image
set(handles.axes,'HandleVisibility','on');
axes(handles.axes);
hold off;
% Display the grayscale image
if isfield(handles,image)
    set(handles.image,'CData',imread([handles.pathname handles.filename]));
else
    handles.image=imshow(imread([handles.pathname handles.filename]));
    hold on
end
% Place the image at the bottom of the stack
uistack(handles.image,'bottom');

set(handles.FileText,'String',[handles.pathname handles.filename]);
files = dir([handles.pathname '*.' handles.filename(end-2:end)]);
[junk, index] = sort(datenum(char(files(:).date)));
set(handles.files_list,'String',{files(index).name});

% Remove the extension of the file (allowing to have dots in the name)
handles.fileroot = handles.filename(1:max(findstr(handles.filename,'.'))-1);

guidata(hObject, handles);  % Update handles structure

% Test existence of the data file
if (exist([handles.pathname handles.fileroot '.mat'],'file') == 2)
    opendata(hObject,eventdata,[handles.pathname handles.fileroot '.mat']);
end

function opendata(hObject,eventdata,file)

handles = guidata(hObject); % Retrieve handles structure

load(file); % Load saved data file

if ~exist('tissue','var'), return; end;
handles.tissue = tissue; % Load tissue
handles.tissueplot  = replot(handles.tissue,handles.tissueplot,...
    'o-','MarkerEdgeColor',[1 1 1],'Linewidth',handles.defaultlinewidth,...
    'MarkerFaceColor',[1 1 1],'Markersize',handles.defaultmarkersize);  % Replot tissue

% Attribute events to graphic objects
set([handles.tissueplot.vertices],'ButtonDownFcn',@select);
set([handles.tissueplot.edges],'ButtonDownFcn',@select);
set([handles.tissueplot.cells],'ButtonDownFcn',@select);
clear tissue

guidata(hObject, handles); % Update handles structure
update_indicators(hObject); %Update indicators

% --------------------------------------------------------------------
function SaveData_Callback(hObject, eventdata, handles)

save_file = [handles.pathname handles.fileroot '.mat'];
tissue = handles.tissue;

save(save_file,'tissue'); % Save tissue object    
clear tissue;

% --------------------------------------------------------------------
function Quit_Callback(hObject, eventdata, handles)

if strcmp(questdlg('Really quit?','Quit Cell Analysis software'),'Yes')
  delete([findobj(0,'tag','sc_importexport') gcbf])
end

% --------------------------------------------------------------------
function About_Callback(hObject, eventdata, handles)

helpdlg(char('Tissue Analysis Software','Version 0.1.',...
    'Last modification: 2009-10-20'),'About');

%--------------------------------------------------------------------
% --- Executes on selection change in files_list.
function files_list_Callback(hObject, eventdata, handles)

% Read the image
string = get(hObject,'String');
im = imread([handles.pathname string{get(hObject,'Value')}]);
set(handles.axes2,'HandleVisibility','on');
axes(handles.axes2);
imshow(im);
axes(handles.axes);

% --- Executes on button press in nextimage.
function nextimage_Callback(hObject, eventdata, handles)

string = get(handles.files_list,'String');
handles.filename = string{get(handles.files_list,'Value')};
guidata(hObject,handles); % Update handles structure
openimage(hObject,eventdata); 
%--------------------------------------------------------------------
function create_Callback(hObject, eventdata, handles)

if get(handles.create,'Value')
    zoom off
    pan off
    % Attribute the functions to the click events
    set([handles.tissueplot.vertices],'ButtonDownFcn',@select);
    set([handles.tissueplot.edges],'ButtonDownFcn','');
    set([handles.tissueplot.cells],'HitTest','off');
    set(handles.image,'ButtonDownFcn',@newVertex)
    set(handles.buttons,'Enable','off');
    set(handles.create,'Enable','on');
else
    select(hObject, eventdata)
    % Attribute the functions to the click events
    set([handles.tissueplot.vertices],'ButtonDownFcn',@select);
    set([handles.tissueplot.edges],'ButtonDownFcn',@select);
    set([handles.tissueplot.cells],'HitTest','on');
    set(handles.image,'ButtonDownFcn','')
    set(handles.buttons,'Enable','on');
end

% --- Executes on button press in edit.
function edit_Callback(hObject, eventdata, handles)

if get(handles.edit,'Value')
    functionname  = @move_Start;
    set(handles.buttons,'Enable','off');
    set(handles.edit,'Enable','on');
else
    functionname = @select;
    set(handles.buttons,'Enable','on');
end
set([handles.tissueplot.vertices],'ButtonDownFcn',functionname);
set([handles.tissueplot.edges],'ButtonDownFcn',functionname);

% --- Executes on button press in edit.
function defineCell_Callback(hObject, eventdata, handles)

if get(handles.defineCell,'Value')
    set([handles.tissueplot.edges],'ButtonDownFcn','');
    set([handles.tissueplot.cells],'ButtonDownFcn','');
    set(handles.buttons,'Enable','off');
    set(handles.defineCell,'Enable','on');
else
    set(handles.buttons,'Enable','on');
    set([handles.tissueplot.vertices],'ButtonDownFcn',@select);
    set([handles.tissueplot.edges],'ButtonDownFcn',@select);
    set([handles.tissueplot.cells],'ButtonDownFcn',@select);
end

function autodefineCell_Callback(hObject, eventdata, handles)

if get(handles.autodefineCell,'Value')
    set([handles.tissueplot.vertices],'ButtonDownFcn','');
    set([handles.tissueplot.edges],'ButtonDownFcn','');
    set([handles.tissueplot.cells],'ButtonDownFcn','');
    set(handles.image,'ButtonDownFcn',@autocelledition)
    set(handles.buttons,'Enable','off');
    set(handles.autodefineCell,'Enable','on');
else
    set(handles.buttons,'Enable','on');
    set([handles.tissueplot.vertices],'ButtonDownFcn',@select);
    set([handles.tissueplot.edges],'ButtonDownFcn',@select);
    set([handles.tissueplot.cells],'ButtonDownFcn',@select);
    set(handles.image,'ButtonDownFcn','')
    select(hObject,eventdata)    
end

%%%%%%%%%%%%%%%%%%% Graphical events %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function select(hObject, eventdata)

handles = guidata(hObject); % Retrieve handles structure
nc = 0;
ne = 0;
nv = 0;

if ~isempty(find(handles.tissueplot.vertices==hObject,1)) % If vertex is selected by clicking on it
    nv = find(handles.tissueplot.vertices==hObject,1);
elseif hObject == handles.SelectedVertex  % If vertex is selected by selecting it on the list
    nv = get(hObject,'Value')-1;    
elseif ~isempty(find(handles.tissueplot.edges==hObject,1)) % If edge is selected by clicking on it
    ne = find(handles.tissueplot.edges==hObject,1);
elseif hObject == handles.SelectedEdge % If edge is selected by selecting it on the list
    ne =get(hObject,'Value')-1;    
elseif ~isempty(find(handles.tissueplot.cells==hObject,1)) % If cell is selected by clicking on it
    nc = find(handles.tissueplot.cells==hObject,1);
elseif hObject == handles.SelectedCell  % If cell is selected by selecting it on the list
    nc =get(hObject,'Value')-1;    
end

% If an element is reselected, cancel the selection
if handles.selectedvertex == nv, nv=0; end;
if handles.selectededge == ne, ne=0; end;
if handles.selectedcell == nc, nc=0; end;

% Case of the creation mode with a previously selected vertex
if get(handles.create,'Value') && nv ~=0 && handles.selectedvertex ~= 0
    age = get(handles.age,'Value');
    handles.tissue = addEdge(handles.tissue,[handles.selectedvertex nv 0 age]);
    handles.tissueplot = replot(handles.tissue,handles.tissueplot,'o-',...
        'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',[1 1 1]); % Replot tissue
end

% Update the selected vertices,edges and cells
handles.selectedvertex = nv;
handles.selectededge = ne;
handles.selectedcell = nc;

% Update marker and lines properties
set([handles.tissueplot.vertices],'MarkerSize',handles.defaultmarkersize);
set([handles.tissueplot.edges],'Linewidth',handles.defaultlinewidth);
set([handles.tissueplot.cells],'FaceColor','none');

set(handles.tissueplot.vertices(nonzeros(nv)),'MarkerSize',handles.selectedmarkersize);
set(handles.tissueplot.edges(nonzeros(ne)),'Linewidth',handles.selectedlinewidth);
set(handles.tissueplot.cells(nonzeros(nc)),'FaceColor','r');
uistack([handles.tissueplot.edges],'top');
uistack([handles.tissueplot.vertices],'top');

guidata(hObject, handles); % Update handles structure
update_indicators(hObject); % Update indicators

% Case of the edition mode 
if get(handles.defineCell,'Value')
    celledition(hObject,eventdata);
end

function celledition(hObject,eventdata)

handles = guidata(hObject); % Retrieve handles structure

if handles.selectedvertex~= 0 % Create new cell
    set([handles.tissueplot.vertices],'ButtonDownFcn','');
    handles.firstcellvertex =handles.selectedvertex;
    lastpoint = handles.selectedvertex;
elseif handles.selectededge ~=0
    if ~isempty(handles.newcell)
    lastedge  = handles.newcell(end);
        if lastedge>0
            lastvertex = handles.tissue.e(lastedge,2);
        else
            lastvertex = handles.tissue.e(-lastedge,1);
        end
    else
        lastvertex = handles.firstcellvertex;
    end

    if handles.tissue.e(handles.selectededge,1) == lastvertex
        handles.newcell(end+1) = handles.selectededge;
        lastpoint = handles.tissue.e(handles.selectededge,2);
    else
        handles.newcell(end+1) = -handles.selectededge;
        lastpoint = handles.tissue.e(handles.selectededge,1);
    end
else
    return
end

if ~isempty(handles.newcell) && handles.firstcellvertex == lastpoint
    handles.tissue = addCell(handles.tissue,handles.newcell);
    handles.newcell = [];
    celledge = [];
    for j =1:length(handles.tissue.c{end})
        celledge = [celledge; edge(handles.tissue,handles.tissue.c{end}(j))];    
    end

    handles.tissueplot.cells(end+1)= fill(celledge(:,1),celledge(:,2),'k','EdgeColor','none','FaceColor','r');

    set([handles.tissueplot.edges],'ButtonDownFcn','','LineStyle','-',...
        'Linewidth',handles.defaultlinewidth);
    set([handles.tissueplot.vertices],'ButtonDownFcn',@select,...
        'MarkerSize',handles.defaultmarkersize);
    handles.selectededge =0;
    handles.selectedcell =length(handles.tissueplot.cells);
    guidata(hObject, handles); % Update handles structure
    update_indicators(hObject); % Update indicators
    return
end

iv = handles.tissue.e(abs(handles.newcell),1:2);
ie = [find(handles.tissue.e(:,1) == lastpoint);...
    find(handles.tissue.e(:,2) == lastpoint)];
ne = ie(ie~=handles.selectededge);
set([handles.tissueplot.edges],'ButtonDownFcn','','LineStyle','-');
set(handles.tissueplot.edges(abs(handles.newcell)),'ButtonDownFcn','',...
    'Linewidth',handles.selectedlinewidth);    
set(handles.tissueplot.edges(ne),'ButtonDownFcn',@select,'LineStyle','--');        
set(handles.tissueplot.vertices(iv(:)),'MarkerSize',handles.selectedmarkersize);
guidata(hObject, handles); % Update handles structure

function autocelledition(hObject,eventdata)

handles = guidata(hObject); % Retrieve handles structure

P = get(handles.axes,'CurrentPoint'); % Get position of the pointer
P = P(1,1:2);
dV =handles.tissue.v-repmat(P,size(handles.tissue.v,1),1);
L=dV(:,1).^2+dV(:,2).^2;
[junk firstpoint] = min(L);

V1=handles.tissue.v(firstpoint,:);
ie1 = find(handles.tissue.e(:,1) == firstpoint);
ie2 = find(handles.tissue.e(:,2) == firstpoint);

theta1=zeros(length(ie1),1);
theta2=zeros(length(ie2),1);
for i =1:length(ie1)
    E1 = handles.tissue.vertices(handles.tissue.e(ie1(i),2),:);
    D = dot(V1-P,E1-P);
    C=cross([V1-P 0],[E1-P 0]);
    theta1(i) = sign(C(3))*acos(D/(norm(V1-P)*norm(E1-P)));
end
for i =1:length(ie2)
    E1 = handles.tissue.vertices(handles.tissue.e(ie2(i),1),:);
    D = dot(V1-P,E1-P);
    C=cross([V1-P 0],[E1-P 0]);
    theta2(i) = sign(C(3))*acos(D/(norm(V1-P)*norm(E1-P)));
end

if isempty(theta1) && isempty(theta2), error('No edge found'); end
if isempty(theta2), mode =1;
elseif isempty(theta1), mode =0;
else mode = min(theta1)<min(theta2);
end
    
if mode
    [junk imin] = min(theta1);   
    lastedge = ie1(imin);
    newcell = lastedge;
    lastlastpoint=handles.tissue.e(lastedge,1);
    lastpoint=handles.tissue.e(lastedge,2);
else
    [junk imin] = min(theta2);
    lastedge = ie2(imin);
    newcell = - lastedge;
    lastlastpoint=handles.tissue.e(lastedge,2);
    lastpoint=handles.tissue.e(lastedge,1);
end


while lastpoint~=firstpoint
    %set(handles.tissueplot.edges(lastedge),'Linewidth',2);
    %waitforbuttonpress
    ie1 = find(handles.tissue.e(:,1) == lastpoint);
    ie1 =ie1(ie1~=lastedge);
    ie2 = find(handles.tissue.e(:,2) == lastpoint);
    ie2 = ie2(ie2~=lastedge);

    oldV=handles.tissue.v(lastlastpoint,:);
    V=handles.tissue.v(lastpoint,:);
    v1= (oldV-V)/norm(V-oldV);
    theta0 =atan2(v1(2),v1(1));
    theta1=zeros(length(ie1),1);
    theta2=zeros(length(ie2),1);
    for i =1:length(ie1)
        newV = handles.tissue.vertices(handles.tissue.e(ie1(i),2),:);
        v2= (newV-V)/norm(newV-V);
        theta1(i) = mod(atan2(v2(2),v2(1))-theta0,2*pi);
    end
    
    for i =1:length(ie2)
        newV = handles.tissue.vertices(handles.tissue.e(ie2(i),1),:);
        v2= (newV-V)/norm(newV-V);
        theta2(i) = mod(atan2(v2(2),v2(1))-theta0,2*pi);
    end

    if isempty(theta1) && isempty(theta2), error('No edge found'); end
    if isempty(theta2), mode =1;
    elseif isempty(theta1), mode =0;
    else mode = min(theta1)<min(theta2);
    end

    if mode
        [junk imin] = min(theta1);
        lastedge = ie1(imin);
        newcell = [newcell lastedge];
        lastlastpoint=handles.tissue.e(lastedge,1);
        lastpoint=handles.tissue.e(lastedge,2);
    else
        [junk imin] = min(theta2);
        lastedge = ie2(imin);
        newcell = [newcell -lastedge];
        lastlastpoint=handles.tissue.e(lastedge,2);
        lastpoint=handles.tissue.e(lastedge,1);
        end
    if isempty(lastpoint), error('emptynewpoint'); end
end


handles.tissue = addCell(handles.tissue,newcell);
if cellArea(handles.tissue,length(handles.tissue.c))<0
    handles.tissue.c{length(handles.tissue.c)}=-handles.tissue.c{length(handles.tissue.c)}(end:-1:1);
end

celledge = [];
for j =1:length(handles.tissue.c{end})
    celledge = [celledge; edge(handles.tissue,handles.tissue.c{end}(j))];    
end

set([handles.tissueplot.cells],'FaceColor','none');
handles.tissueplot.cells(end+1)= fill(celledge(:,1),celledge(:,2),'k','EdgeColor','none','FaceColor','r');

handles.selectedcell =length(handles.tissueplot.cells);
guidata(hObject, handles); % Update handles structure
update_indicators(hObject); % Update indicators

function newVertex(hObject, eventdata)

handles = guidata(hObject); % Retrieve handles.structure

P = get(handles.axes,'CurrentPoint'); % Get position of the pointer
handles.tissue = addVertex(handles.tissue,P(1,1:2)); % Add new vertex to tissue
set([handles.tissueplot.vertices],'MarkerSize',handles.defaultmarkersize);
handles.tissueplot.vertices(end+1) = plot(handles.tissue.v(end,1), handles.tissue.v(end,2),'o',...
        'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',[1 1 1],...
        'ButtonDownFcn',@select,'MarkerSize',handles.selectedmarkersize,'LineStyle','none');
    
nv = size(handles.tissue.v,1);
if handles.selectedvertex % If a previous point has been selected
    age = get(handles.age,'Value');
    handles.tissue =addEdge(handles.tissue,[handles.selectedvertex nv 0 age]); % Add new edge to tissue
    [newedge edgecolor] = edge(handles.tissue,size(handles.tissue.e,1));
    handles.tissueplot.edges(end+1) = plot(newedge(:,1), newedge(:,2),...
        'Linewidth',handles.defaultlinewidth,'Color',edgecolor,'Marker','none');    
end
uistack([handles.tissueplot.edges],'top');
uistack([handles.tissueplot.vertices],'top');

handles.selectedvertex = nv;

guidata(hObject, handles); % Update handles structure
update_indicators(hObject, eventdata) % Update indicators

function move_Start(hObject, eventdata)

handles = guidata(hObject); % Retrieve handles structure

if ~isempty(find(handles.tissueplot.vertices==hObject,1)) % If vertex is selected
    handles.selectedvertex = find(handles.tissueplot.vertices==hObject,1);
elseif ~isempty(find(handles.tissueplot.edges==hObject,1))% If edge is selected
    handles.selectededge = find(handles.tissueplot.edges==hObject,1);
end

set([handles.tissueplot.edges],'EraseMode','xor');
set([handles.tissueplot.vertices],'EraseMode','xor');
guidata(hObject, handles); % Update handles structure

% Update motion properties
set(handles.figure,'WindowButtonMotion',@move_Motion,'WindowButtonUp',@move_Finish)

function move_Motion(hObject, eventdata)

handles = guidata(hObject); % Retrieve handles structure
handles = editNetwork(handles);
guidata(hObject, handles); % Update handles structure

function move_Finish(hObject, eventdata)

handles = guidata(hObject); % Retrieve handles structure
handles = editNetwork(handles);
% Remove the motion events
set(handles.figure,'WindowButtonMotion','','WindowButtonUp','')
set([handles.tissueplot.vertices],'EraseMode','normal');
set([handles.tissueplot.edges],'EraseMode','normal');

% Reinitialize constants
handles.selectedvertex = 0;
handles.selectededge = 0;
handles.selectedcell = 0;

guidata(hObject, handles); % Update handles structure

function handles = editNetwork(handles)

P = get(handles.axes,'CurrentPoint'); % Get position of the pointer
axlim=axis;
newP(1,1) = min(max(P(1,1),axlim(1)), axlim(2));
newP(1,2) = min(max(P(1,2),axlim(3)), axlim(4));
nv =[];
ne =[];
nc = [];
if handles.selectedvertex ~=0
    % Retrieve and update vertices list
    vertices = get(handles.tissue,'Vertices');
    vertices(handles.selectedvertex,:) = newP;
    handles.tissue = set(handles.tissue,'Vertices',vertices);

    nv = handles.selectedvertex;
    ne = [find(handles.tissue.e(:,1) ==handles.selectedvertex); ...
        find(handles.tissue.e(:,2) ==handles.selectedvertex)];
    nc = neighbors(handles.tissue,ne);
else
    % Retrieve and update edge list
    vertices = get(handles.tissue,'Vertices');
    edges = get(handles.tissue,'Edges');
    V1 = vertices(edges(handles.selectededge,1),:);
    V2 = vertices(edges(handles.selectededge,2),:);
    v1 = (V1-newP)/norm(V1-newP);
    v2 = (V2-newP)/norm(V2-newP);
    vperp = [-v1(2) v1(1)];
    alpha = atan2(dot(v2,vperp),dot(v2,v1));
    if alpha>0
        theta = alpha-pi;
    else
        theta = pi+alpha;
    end
    edges(handles.selectededge,3)=theta;
    handles.tissue = set(handles.tissue,'Edges',edges);
    ne = handles.selectededge;
    nc = neighbors(handles.tissue,ne);
end
%guidata(hObject, handles); % Update handles structure

for i=nv
    set(handles.tissueplot.vertices(i),'XData',handles.tissue.v(i,1),'YData',handles.tissue.v(i,2));
end
for i=1:length(ne)
    E =edge(handles.tissue,ne(i));
    set(handles.tissueplot.edges(ne(i)),'XData',E(:,1),'YData',E(:,2));
end
for i=1:length(nc)
    celledge = [];
    for j =1:length(handles.tissue.c{nc(i)})
        celledge = [celledge; edge(handles.tissue,handles.tissue.c{nc(i)}(j))];    
    end
    set(handles.tissueplot.cells(nc(i)),'XData',celledge(:,1),'YData',celledge(:,2));
end


% --- Executes on button press in clear.
function clear_Callback(hObject, eventdata, handles)

if handles.selectedvertex
    handles.tissue=removeVertex(handles.tissue, handles.selectedvertex);
end
if handles.selectededge
    handles.tissue=removeEdge(handles.tissue, handles.selectededge);
end
if handles.selectedcell
    handles.tissue=removeCell(handles.tissue, handles.selectedcell);
end
% Replot cell
handles.tissueplot = replot(handles.tissue,handles.tissueplot,...
    'o-','MarkerEdgeColor',[1 1 1],'Linewidth',handles.defaultlinewidth,...
    'MarkerFaceColor',[1 1 1],'Markersize',handles.defaultmarkersize);
handles.selectedvertex = 0;
handles.selectededge = 0;
handles.selectedcell = 0;

guidata(hObject, handles); % Update handles structure
update_indicators(hObject); % Update indicators

function update_indicators(hObject, eventdata)

handles = guidata(hObject); % Retrieve handles structure

nv = size(handles.tissue.v,1);
vertexstring(1) = {'No Vertex'};
vertexstring(2:nv+1)=cellstr([repmat('Vertex ',nv,1) num2str((1:nv)')]);
set(handles.SelectedVertex,'Value',handles.selectedvertex+1,'String',vertexstring);

ne = size(handles.tissue.e,1);
edgestring(1) = {'No Edge'};
edgestring(2:ne+1)=cellstr([repmat('Edge ',ne,1) num2str((1:ne)')]);
set(handles.SelectedEdge,'Value',handles.selectededge+1,'String',edgestring);

nc = length(handles.tissue.c);
cellstring(1) = {'No Cell'};
cellstring(2:nc+1)=cellstr([repmat('Cell ',nc,1) num2str((1:nc)')]);
set(handles.SelectedCell,'Value',handles.selectedcell+1,'String',cellstring);

if handles.selectededge ~= 0
    edges = get(handles.tissue,'Edges');
    set(handles.age,'Value',edges(handles.selectededge,4));
    set(handles.agelabel,'String',num2str(edges(handles.selectededge,4)),...
        'BackgroundColor',coloredge(edges(handles.selectededge,4)));
end

% --- Executes on slider movement.
function age_Callback(hObject, eventdata, handles)

age = round(get(hObject,'Value'));
set(handles.agelabel,'String',num2str(age),'BackgroundColor',coloredge(age));

if handles.selectededge ~=0
    edges = get(handles.tissue,'Edges');	% Retrieve list of edges
    edges(handles.selectededge,4) = age;    % Modify age of the selected edge
    handles.tissue = set(handles.tissue,'Edges',edges); % Modify tissue object
    guidata(hObject,handles);	% Update handles structure
    set(handles.tissueplot.edges(handles.selectededge),'Color',coloredge(age));
end
