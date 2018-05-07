function varargout = MONOVIEW(varargin)
% Monoview is a simple interface written in Matlab using its graphical user 
% interface (GUI) module, GUIDE. The tool enables the user to click on one 
% main image (Cam1) and to be immediately directed towards other photos 
% where the clicked point is visible (Cam2). Furthermore, by using the 
% monoplotting concept, the user will also be able to inquire the 3D 
% coordinates of the selected point. Other features of the tool include a 
% distance calculator between two points in one image, and a vectorizer 
% which can be used to digitize polygonal features from an image.
% 
% Tape "guide" on Matlab command window and open MONOVIEW.fig to see 
% the graphical interface.
%
% Written by Arnadi MURTIYOSO
% ifp Stuttgart & INSA de Strasbourg
% last update: 07/07/15
%

format long g
clc
% MONOVIEW MATLAB code for MONOVIEW.fig
%      MONOVIEW, by itself, creates a new MONOVIEW or raises the existing
%      singleton*.
%
%      H = MONOVIEW returns the handle to a new MONOVIEW or the handle to
%      the existing singleton*.
%
%      MONOVIEW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MONOVIEW.M with the given input arguments.
%
%      MONOVIEW('Property','Value',...) creates a new MONOVIEW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MONOVIEW_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MONOVIEW_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MONOVIEW

% Last Modified by GUIDE v2.5 30-Jun-2015 10:18:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MONOVIEW_OpeningFcn, ...
                   'gui_OutputFcn',  @MONOVIEW_OutputFcn, ...
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
% End initialization code - DO NOT EDIT


% --- Executes just before MONOVIEW is made visible.
function MONOVIEW_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MONOVIEW (see VARARGIN)

% Choose default command line output for MONOVIEW
handles.output = hObject;
set(handles.PointSelect,'Enable', 'Off');
set(handles.PointSelect2,'Enable', 'Off');
set(handles.zoomto,'Enable', 'Off');
E2={};
Dxf={};
A={};
handles.A = A;
handles.E2 = E2;
handles.Dxf=Dxf;
%set(handles.zoomto,'Enable', 'Off');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MONOVIEW wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MONOVIEW_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% =========================================================================
%Browse directory button
% =========================================================================
% --- Executes on button press in BrwsList1.
function BrwsList1_Callback(hObject, eventdata, handles)
% hObject    handle to BrwsList1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set (handles.PointSelect, 'Enable', 'Off');
handles.output = hObject;
    [fn] = uigetdir('', 'Select directory');
complete = strcat (fn);
handles.DirFile = complete;
guidata(hObject, handles);

filepattern = fullfile(complete,'/*.*');
ImageFiles = dir(filepattern);

ListofImageNames = {};
for Index = 1:length(ImageFiles)
    baseFileName = ImageFiles(Index).name;
    [folder, name, extension] = fileparts (baseFileName);
    extension = upper(extension);
   
    if strcmp (extension, '.JPG') == true %&& isempty(regexp(name,'Cam1')) == false
        ListofImageNames = [ListofImageNames baseFileName];
    else
    end
end
set(handles.listCam1, 'string', ListofImageNames);

index_selected = get(handles.listCam1,'Value');
file_list = get(handles.listCam1,'String');
    filename = file_list{index_selected};
fulladdress = strcat(complete,'/',filename);

[~, file, ~] = fileparts (filename);
file1 = file;
save('file1.mat', 'file1');

J = imread(fulladdress);
axes(handles.Cam1);
handles.Cam1 = imshow(J,[]);

% =========================================================================
%Cam 1 peripherals
% =========================================================================
% --- Executes on selection change in listCam1.
function listCam1_Callback(hObject, eventdata, handles)

index_selected = get(handles.listCam1,'Value');
file_list = get(handles.listCam1,'String');
    filename = file_list{index_selected};
% handles.index1 = index_selected;

DirFile=handles.DirFile;
complete = strcat(DirFile,'/',filename);

J = imread(complete);
axes(handles.Cam1);
handles.Cam1 = imshow(J,[]);

[~, file, ~] = fileparts (filename);
file1 = file;
save('file1.mat', 'file1');
% hObject    handle to listCam1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listCam1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listCam1


% --- Executes during object creation, after setting all properties.
function listCam1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listCam1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listCam2.
function listCam2_Callback(hObject, eventdata, handles)

index_selected = get(handles.listCam2,'Value');
file_list = get(handles.listCam2,'String');
    filename = file_list{index_selected};
% handles.index1 = index_selected;

DirFile=handles.DirFile;
complete = strcat(DirFile,'/',filename);

J = imread(complete);
axes(handles.Cam2);
handles.Cam2 = imshow(J,[]);

[~, file, ~] = fileparts (filename);
file2 = file;
save('file2.mat', 'file2');
% hObject    handle to listCam2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listCam2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listCam2


% --- Executes during object creation, after setting all properties.
function listCam2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listCam2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editEO1_Callback(hObject, eventdata, handles)
% hObject    handle to editEO1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editEO1 as text
%        str2double(get(hObject,'String')) returns contents of editEO1 as a double


% --- Executes during object creation, after setting all properties.
function editEO1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editEO1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in BrwsEO1.
function BrwsEO1_Callback(hObject, eventdata, handles)
handles.output = hObject;
[fn pn] = uigetfile('*.txt', 'Select EO file(*.txt)');
complete = strcat (pn,fn);
set (handles.editEO1, 'string', complete);

FID = fopen (complete);
datacell = textscan (FID, '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s', 'HeaderLines', 2, 'CollectOutput', 1);
fclose(FID);
A = datacell{1};
save('A.mat','A');

guidata(handles.BrwsEO1, handles);
% hObject    handle to BrwsEO1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function editIO1_Callback(hObject, eventdata, handles)
% hObject    handle to editIO1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editIO1 as text
%        str2double(get(hObject,'String')) returns contents of editIO1 as a double


% --- Executes during object creation, after setting all properties.
function editIO1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editIO1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in BrwsIO1.
function BrwsIO1_Callback(hObject, eventdata, handles)
handles.output = hObject;
[fn pn] = uigetfile('*.xml', 'Select IO file(*.xml)');
complete = strcat (pn,fn);
set (handles.editIO1, 'string', complete);
pix = str2num(get(handles.pixel,'String'));
 xdoc = xmlread(complete);
 xmlwrite(xdoc);
 
 import javax.xml.xpath.*
 factory = XPathFactory.newInstance;
 xpath = factory.newXPath;
 
 expression = xpath.compile('calibration/width');
 pixwNode = expression.evaluate(xdoc, XPathConstants.NODE);
 pixw = pixwNode.getTextContent;
 pixw = str2num (pixw);
 
 expression = xpath.compile('calibration/height');
 pixhNode = expression.evaluate(xdoc, XPathConstants.NODE);
 pixh = pixhNode.getTextContent;
 pixh = str2num (pixh);
 
 
 expression = xpath.compile('calibration/cx');
 cxNode = expression.evaluate(xdoc, XPathConstants.NODE);
 cx = cxNode.getTextContent;
 cx = str2num (cx);
 
 expression = xpath.compile('calibration/cy');
 cyNode = expression.evaluate(xdoc, XPathConstants.NODE);
 cy = cyNode.getTextContent;
 cy = str2num (cy);
 
 x01=cx*pix;
 x01=x01-((pix*pixw)/2);
 y01=cy*pix;
 y01=-y01+((pix*pixh)/2);
 
 expression = xpath.compile('calibration/fx');
 fxNode = expression.evaluate(xdoc, XPathConstants.NODE);
 fx = fxNode.getTextContent;
 fx = str2num (fx);
 
 expression = xpath.compile('calibration/fy');
 fyNode = expression.evaluate(xdoc, XPathConstants.NODE);
 fy = fyNode.getTextContent;
 fy = str2num (fy);
 
 fx = fx*pix;
 fy = fy*pix;
 
 c1 = (fx+fy)/2;
 
 
 expression = xpath.compile('calibration/k1');
 K1Node = expression.evaluate(xdoc, XPathConstants.NODE);
 K1 = K1Node.getTextContent;
 K1 = str2num (K1);
 
 K1 = K1/fx^2;
 
 expression = xpath.compile('calibration/k2');
 K2Node = expression.evaluate(xdoc, XPathConstants.NODE);
 K2 = K2Node.getTextContent;
 K2 = str2num (K2);
 
 K2 = K2/fx^4;
 
 expression = xpath.compile('calibration/k3');
 K3Node = expression.evaluate(xdoc, XPathConstants.NODE);
 K3 = K3Node.getTextContent;
 K3 = str2num (K3);
 
 K3 = K3/fx^6;
 
 expression = xpath.compile('calibration/p1');
 P1Node = expression.evaluate(xdoc, XPathConstants.NODE);
 P1 = P1Node.getTextContent;
 P1 = str2num (P1);
 
 expression = xpath.compile('calibration/p2');
 P2Node = expression.evaluate(xdoc, XPathConstants.NODE);
 P2 = P2Node.getTextContent;
 P2 = str2num (P2);
 

%  
% pixw = str2num(get(handles.pixwidth,'String'));
% pixh = str2num(get(handles.pixheight,'String'));
% pix = str2num(get(handles.pixel,'String'));
%  
%  x01 = x01*pix;
%  y01 = y01*pix;
%  c1 = c1*pix;
%  
%  x01 = x01-((pix*pixw)/2);
%  y01 = -y01+((pix*pixh)/2);
    
 %  RNode = expression.evaluate(xdoc, XPathConstants.NODE);
%  R = RNode.getTextContent;
%  R = str2num (R)
 
 handles.pixw = pixw;
 handles.pixh = pixh;
 handles.x01 = x01;
 handles.y01 = y01;
 handles.c1 = c1;
 handles.K1 = K1;
 handles.K2 = K2;
 handles.K3 = K3;
 handles.P1 = P1;
 handles.P2 = P2;
 
guidata(handles.BrwsIO1, handles);
% hObject    handle to BrwsIO1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function pixel_Callback(hObject, eventdata, handles)
% hObject    handle to pixel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pixel as text
%        str2double(get(hObject,'String')) returns contents of pixel as a double


% --- Executes during object creation, after setting all properties.
function pixel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pixel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pixwidth_Callback(hObject, eventdata, handles)
% hObject    handle to pixwidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pixwidth as text
%        str2double(get(hObject,'String')) returns contents of pixwidth as a double


% --- Executes during object creation, after setting all properties.
function pixwidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pixwidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pixheight_Callback(hObject, eventdata, handles)
% hObject    handle to pixheight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pixheight as text
%        str2double(get(hObject,'String')) returns contents of pixheight as a double


% --- Executes during object creation, after setting all properties.
function pixheight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pixheight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PointSelect.
function PointSelect_Callback(hObject, eventdata, handles)
set(handles.zoomto,'Enable', 'On');
set(handles.PointSelect2, 'Enable', 'On');
[x1,y1] = ginput (1);
handles.X = x1;
handles.Y = y1;
hold on;
plot(x1,y1,'y+');
hold on;
guidata(handles.PointSelect,handles);
pixw = handles.pixw;
pixh = handles.pixh;
% pixw = str2num(get(handles.pixwidth,'String'));
% pixh = str2num(get(handles.pixheight,'String'));
pix = str2num(get(handles.pixel,'String'));

x1=x1*pix;
x1=x1-((pix*pixw)/2);
y1=y1*pix;
y1=-y1+((pix*pixh)/2);

storedStructure = load('Cam1.mat');
Cam1 = storedStructure.Cam1;

% hObject    handle to PointSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%interior orientation
c1 = handles.c1;
x01 = handles.x01;
y01 = handles.y01;
K1 = handles.K1;
K2 = handles.K2;
K3 = handles.K3;
P1 = handles.P1;
P2 = handles.P2;


%exterior orientation
storedStructure2 = load('file1.mat');
file1 = storedStructure2.file1;
file1 = strcat(file1,'.JPG');
 row = Cam1(strcmp(file1,Cam1(:,1)),:);
% row = find (ismember(Cam1,file1),1)

X01 = cell2mat(row (1,2));
Y01 = cell2mat(row (1,3));
Z01 = cell2mat(row (1,4));
omega1 = cell2mat(row (1,5));
phi1 = cell2mat(row (1,6));
kappa1 = cell2mat(row (1,7));

%rotation matrix
r1=rotation2(omega1,phi1,kappa1);

%bilinear interpolation for ground altitude 
%(1st initial value: X Y of camera)
% Z = 250;
storedStructure2 = load('DSM.mat');
DSM = storedStructure2.DSM;

X0tiff = handles.X0tiff;
Y0tiff = handles.Y0tiff;
ResTiff = handles.ResTiff;

Z = bilineinter (DSM,X0tiff,Y0tiff,ResTiff,X01,Y01);

x = x1-x01;
y = y1-y01;

%calibration of image coordinates
%radial corrections
r = sqrt(x^2+y^2);
% dr = K1*r^3-K2*r^5+K3*r^7;
% drx = x1*(dr/r);
% dry = y1*(dr/r);
drx = x*(K1*r^2+K2*r^4+K3*r^6);
dry = y*(K1*r^2+K2*r^4+K3*r^6);

%decentering corrections
ddx = P1*(3*x^2+y^2)+2*P2*x*y;
ddy = P2*(x^2+3*y^2)+2*P1*x*y;

x = x - drx - ddx;
y = y - dry - ddy;

%loop to iterate Z
for i=1:3
    %collinearity
    Xx = (Z-Z01)*((r1(1,1)*x + r1(2,1)*y - r1(3,1)*c1)/(r1(1,3)*x + r1(2,3)*y - r1(3,3)*c1));
    Xx = sprintf('%.4f',Xx);
    Xx = str2double(Xx);
    X01 = sprintf('%.4f',X01);
    X01 = str2double(X01);
    X = X01 + Xx;
    
    Yy = (Z-Z01)*((r1(1,2)*x + r1(2,2)*y - r1(3,2)*c1)/(r1(1,3)*x + r1(2,3)*y - r1(3,3)*c1));
    Yy = sprintf('%.4f',Yy);
    Yy = str2double(Yy);
    Y01 = sprintf('%.4f',Y01);
    Y01 = str2double(Y01);
    Y = Y01 + Yy;
    
    Z = bilineinter (DSM,X0tiff,Y0tiff,ResTiff,X,Y);
end

set (handles.editXres, 'string', sprintf('%0.3f',X));
set (handles.editYres, 'string', sprintf('%0.3f',Y));
set (handles.editZres, 'string', sprintf('%0.3f',Z));

%cam 2
%interior orientation
c2 = handles.c1;
x02 = handles.x01;
y02 = handles.y01;

[r,c] = size(Cam1);
for i = 1:r
    %exterior orientation
    X02 = cell2mat(Cam1(i,2));
    Y02 = cell2mat(Cam1(i,3));
    Z02 = cell2mat(Cam1(i,4));
    omega2 = cell2mat(Cam1(i,5));
    phi2 = cell2mat(Cam1(i,6));
    kappa2 = cell2mat(Cam1(i,7));
%rotation matrix
r2 = rotation2(omega2, phi2, kappa2);

%collinearity
Xt2 = X-X02;
Yt2 = Y-Y02;
Zt2 = Z-Z02;

x2 = x02 - c2*((r2(1,1)*Xt2 + r2(1,2)*Yt2 + r2(1,3)*Zt2)/(r2(3,1)*Xt2 + r2(3,2)*Yt2 + r2(3,3)*Zt2)) ;
y2 = y02 - c2*((r2(2,1)*Xt2 + r2(2,2)*Yt2 + r2(2,3)*Zt2)/(r2(3,1)*Xt2 + r2(3,2)*Yt2 + r2(3,3)*Zt2)) ;

%calibration of image coordinates
%radial corrections
r = sqrt(x2^2+y2^2);
% dr = K1*r^3-K2*r^5+K3*r^7;
% drx = x1*(dr/r);
% dry = y1*(dr/r);
drx = x2*(K1*r^2+K2*r^4+K3*r^6);
dry = y2*(K1*r^2+K2*r^4+K3*r^6);

%decentering corrections
ddx = P1*(3*x2^2+y2^2)+2*P2*x2*y2;
ddy = P2*(x2^2+3*y2^2)+2*P1*x2*y2;

x2 = x2 + drx + ddx;
y2 = y2 + dry + ddy;


x2 = x2+((pix*pixw)/2);
x2 = x2/pix;
y2 = -(y2-((pix*pixh)/2));
y2 = y2/pix;

Cam1{i,8}=(x2);
Cam1{i,9}=(y2);

end


save('Cam1.mat', 'Cam1');
ListofImageNames = {};
for Index = 1:r
    baseFileName = Cam1{Index,1};
    x2 = cell2mat(Cam1(Index,8));
    y2 = cell2mat(Cam1(Index,9));
    [folder, name, extension] = fileparts (baseFileName);
    extension = upper(extension);
   
    if y2>0 && x2>0 && y2<pixh && x2<pixw
        ListofImageNames = [ListofImageNames baseFileName];
    else
    end
end
set(handles.listCam2, 'string', ListofImageNames)

DirFile = handles.DirFile;

index_selected = get(handles.listCam2,'Value');
file_list = get(handles.listCam2,'String');
    filename = file_list{index_selected};
fulladdress = strcat(DirFile,'\',filename);

[~, file, ~] = fileparts (filename);
file2 = file;
save('file2.mat', 'file2');

J = imread(fulladdress);
axes(handles.Cam2);
handles.Cam2 = imshow(J,[]);


% hObject    handle to PointSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in setCam.
function setCam_Callback(hObject, eventdata, handles)
set (handles.PointSelect, 'Enable', 'On');

DirFile = handles.DirFile;
ImageFiles = dir(DirFile);

storedStructure = load('A.mat');
A = storedStructure.A;

Cam1 = {};
for Index = 1:length(ImageFiles)
    baseFileName = ImageFiles(Index).name;
    [folder, name, extension] = fileparts (baseFileName);
    extension = upper(extension);
   
    if strcmp (extension, '.JPG') == true %&& isempty(regexp(name,'Cam1')) == false
        Cam1 = [Cam1 strcat(name,extension)];
    else
    end
end
Cam1 = transpose (Cam1);

% populate cell array with EO
for i = 1:length (Cam1)
    name = Cam1(i,1);
    r = find (ismember(A,name),1);

    Cam1{i,2} = str2double(A (r,2));
    Cam1{i,3} = str2double(A (r,3));
    Cam1{i,4} = str2double(A (r,4));
    Cam1{i,5} = str2double(A (r,5));
    Cam1{i,6} = str2double(A (r,6));
    Cam1{i,7} = str2double(A (r,7));
end

%load DSM
dsmpath = handles.dsmpath;
DSM = imread(dsmpath);
worldFileName = getworldfilename(dsmpath);
R = worldfileread(worldFileName, 'planar', size(DSM));
X0tiff = R.XWorldLimits(1,1);
Y0tiff = R.YWorldLimits(1,2);
ResTiff = R.CellExtentInWorldX;

handles.X0tiff = X0tiff;
handles.Y0tiff = Y0tiff;
handles.ResTiff = ResTiff;

save('Cam1.mat','Cam1');
save('DSM.mat','DSM');

guidata(handles.setCam, handles);
% hObject    handle to setCam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in zoomto.
function zoomto_Callback(hObject, eventdata, handles)
storedStructure = load('Cam1.mat');
Cam2 = storedStructure.Cam1;

storedStructure = load('file2.mat');
file2 = storedStructure.file2;
file2 = strcat(file2,'.JPG');

row = Cam2(strcmp(file2,Cam2(:,1)),:);
T = [row(1,8) row(1,9)];
T = cell2mat(T);
x2 = T(1,1);
y2 = T(1,2);

h = handles.Cam2;

set(handles.Cam2, 'xLim', [x2-300 x2+300])
set(handles.Cam2, 'yLim', [y2-200 y2+200])

hold on;
plot(x2,y2,'rx');
hold on;
% hObject    handle to zoomto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function editDSM_Callback(hObject, eventdata, handles)
% hObject    handle to editDSM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDSM as text
%        str2double(get(hObject,'String')) returns contents of editDSM as a double


% --- Executes during object creation, after setting all properties.
function editDSM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDSM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in BrwsDSM.
function BrwsDSM_Callback(hObject, eventdata, handles)
handles.output = hObject;
[fn pn] = uigetfile('*.tif', 'Select DSM file(*.tif)');
complete = strcat (pn,fn);
set (handles.editDSM, 'string', complete);

handles.dsmpath = complete;

guidata(handles.BrwsDSM, handles);

% hObject    handle to BrwsDSM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function editXres_Callback(hObject, eventdata, handles)
% hObject    handle to editXres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editXres as text
%        str2double(get(hObject,'String')) returns contents of editXres as a double


% --- Executes during object creation, after setting all properties.
function editXres_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editXres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editYres_Callback(hObject, eventdata, handles)
% hObject    handle to editYres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editYres as text
%        str2double(get(hObject,'String')) returns contents of editYres as a double


% --- Executes during object creation, after setting all properties.
function editYres_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editYres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editZres_Callback(hObject, eventdata, handles)
% hObject    handle to editZres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editZres as text
%        str2double(get(hObject,'String')) returns contents of editZres as a double


% --- Executes during object creation, after setting all properties.
function editZres_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editZres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Distance_Callback(hObject, eventdata, handles)
% hObject    handle to Distance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Distance as text
%        str2double(get(hObject,'String')) returns contents of Distance as a double


% --- Executes during object creation, after setting all properties.
function Distance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Distance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupDist.
function popupDist_Callback(hObject, eventdata, handles)
% hObject    handle to popupDist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupDist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupDist


% --- Executes during object creation, after setting all properties.
function popupDist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupDist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in toggleDist.
function toggleDist_Callback(hObject, eventdata, handles)
i=1;
A=[];
while get(hObject,'Value')
   [x, y] = ginput(1);
   A(i,1) = x;
   A(i,2) = y;
   i = i+1;
   if i==3
       set (hObject, 'Value',0)
       break        
   end
end

storedStructure = load('Cam1.mat');
Cam1 = storedStructure.Cam1;

cam = get(handles.popupDist,'Value');

c = handles.c1;
x0 = handles.x01;
y0 = handles.y01;
K1 = handles.K1;
K2 = handles.K2;
K3 = handles.K3;
P1 = handles.P1;
P2 = handles.P2;

xp1 = A(1,1);
yp1 = A(1,2);
xp2 = A(2,1);
yp2 = A(2,2);

X = [xp1; xp2];
Y = [yp1; yp2];
if cam == 1
    storedStructure2 = load('file1.mat');
    file1 = storedStructure2.file1;
    file1 = strcat(file1,'.JPG');
    row= Cam1(strcmp(file1,Cam1(:,1)),:);

    X0 = cell2mat(row (1,2));
    Y0 = cell2mat(row (1,3));
    Z0 = cell2mat(row (1,4));
    omega = cell2mat(row (1,5));
    phi = cell2mat(row (1,6));
    kappa = cell2mat(row (1,7));
    
    axes(handles.Cam1);
    hold on
    plot (X,Y,'-+r',...
        'LineWidth',2)
    hold on
elseif cam == 2
    storedStructure2 = load('file2.mat');
    file2 = storedStructure2.file2;
    file2 = strcat(file2,'.JPG');
    row= Cam1(strcmp(file2,Cam1(:,1)),:);

    X0 = cell2mat(row (1,2));
    Y0 = cell2mat(row (1,3));
    Z0 = cell2mat(row (1,4));
    omega = cell2mat(row (1,5));
    phi = cell2mat(row (1,6));
    kappa = cell2mat(row (1,7));
    
    axes(handles.Cam2);
    hold on
    plot (X, Y,'-+r',...
        'LineWidth',2)
    hold on
else
    disp('An error has occured')
end

pixw = handles.pixw;
pixh = handles.pixh;
% pixw = str2num(get(handles.pixwidth,'String'));
% pixh = str2num(get(handles.pixheight,'String'));
pix = str2num(get(handles.pixel,'String'));

xp1 = (xp1*pix)-((pix*pixw)/2);
yp1 = -(yp1*pix)+((pix*pixh)/2);
xp2 = (xp2*pix)-((pix*pixw)/2);
yp2 = -(yp2*pix)+((pix*pixh)/2);


storedStructure2 = load('DSM.mat');
DSM = storedStructure2.DSM;

X0tiff = handles.X0tiff;
Y0tiff = handles.Y0tiff;
ResTiff = handles.ResTiff;
%=======
%POINT 1
%=======
%bilinear interpolation for ground altitude 
%(1st initial value: X Y of camera)
% Z = 250;
Z = bilineinter (DSM,X0tiff,Y0tiff,ResTiff,X0,Y0);

r1=rotation2(omega,phi,kappa);

x = xp1 -x0;
y = yp1 -y0;

%calibration of image coordinates
%radial corrections
r = sqrt(x^2+y^2);
% dr = K1*r^3-K2*r^5+K3*r^7;
% drx = x1*(dr/r);
% dry = y1*(dr/r);
drx = x*(K1*r^2+K2*r^4+K3*r^6);
dry = y*(K1*r^2+K2*r^4+K3*r^6);

%decentering corrections
ddx = P1*(3*x^2+y^2)+2*P2*x*y;
ddy = P2*(x^2+3*y^2)+2*P1*x*y;

x = xp1 -x0 - drx - ddx;
y = yp1 -y0 - dry - ddy;

%Point 1 coordinates
for i=1:3
    %collinearity
    Xx = (Z-Z0)*((r1(1,1)*x + r1(2,1)*y - r1(3,1)*c)/(r1(1,3)*x + r1(2,3)*y - r1(3,3)*c));
    Xx = sprintf('%.4f',Xx);
    Xx = str2double(Xx);
    X0 = sprintf('%.4f',X0);
    X0 = str2double(X0);
    X = X0 + Xx;
    
    Yy = (Z-Z0)*((r1(1,2)*x + r1(2,2)*y - r1(3,2)*c)/(r1(1,3)*x + r1(2,3)*y - r1(3,3)*c));
    Yy = sprintf('%.4f',Yy);
    Yy = str2double(Yy);
    Y0 = sprintf('%.4f',Y0);
    Y0 = str2double(Y0);
    Y = Y0 + Yy;
    
    Z = bilineinter (DSM,X0tiff,Y0tiff,ResTiff,X,Y);
end
Xp1 = X;
Yp1 = Y;
Zp1 = Z;

%=======
%POINT 2
%=======
%bilinear interpolation for ground altitude 
%(1st initial value: X Y of camera)
% Z = 250;
Z = bilineinter (DSM,X0tiff,Y0tiff,ResTiff,X0,Y0);

r1=rotation2(omega,phi,kappa);

x = xp2 -x0;
y = yp2 -y0;
%calibration of image coordinates
%radial corrections
r = sqrt(x^2+y^2);
% dr = K1*r^3-K2*r^5+K3*r^7;
% drx = x1*(dr/r);
% dry = y1*(dr/r);
drx = x*(K1*r^2+K2*r^4+K3*r^6);
dry = y*(K1*r^2+K2*r^4+K3*r^6);

%decentering corrections
ddx = P1*(3*x^2+y^2)+2*P2*x*y;
ddy = P2*(x^2+3*y^2)+2*P1*x*y;

x = xp2 -x0 - drx - ddx;
y = yp2 -y0 - dry - ddy;

%Point 1 coordinates
for i=1:3
    %collinearity
    Xx = (Z-Z0)*((r1(1,1)*x + r1(2,1)*y - r1(3,1)*c)/(r1(1,3)*x + r1(2,3)*y - r1(3,3)*c));
    Xx = sprintf('%.4f',Xx);
    Xx = str2double(Xx);
    X0 = sprintf('%.4f',X0);
    X0 = str2double(X0);
    X = X0 + Xx;
    
    Yy = (Z-Z0)*((r1(1,2)*x + r1(2,2)*y - r1(3,2)*c)/(r1(1,3)*x + r1(2,3)*y - r1(3,3)*c));
    Yy = sprintf('%.4f',Yy);
    Yy = str2double(Yy);
    Y0 = sprintf('%.4f',Y0);
    Y0 = str2double(Y0);
    Y = Y0 + Yy;
    
    Z = bilineinter (DSM,X0tiff,Y0tiff,ResTiff,X,Y);
end
Xp2 = X;
Yp2 = Y;
Zp2 = Z;



d3 = sqrt((Xp1-Xp2)^2+(Yp1-Yp2)^2+(Zp1-Zp2)^2);
d2 = sqrt((Xp1-Xp2)^2+(Yp1-Yp2)^2);
d1 = Zp2-Zp1;

set (handles.Distance, 'string', sprintf('%0.3f',d3));
set (handles.Distance2D, 'string', sprintf('%0.3f',d2));
set (handles.DistanceHeight, 'string', sprintf('%0.3f',d1));
% hObject    handle to toggleDist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of toggleDist



function Distance2D_Callback(hObject, eventdata, handles)
% hObject    handle to Distance2D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Distance2D as text
%        str2double(get(hObject,'String')) returns contents of Distance2D as a double


% --- Executes during object creation, after setting all properties.
function Distance2D_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Distance2D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DistanceHeight_Callback(hObject, eventdata, handles)
% hObject    handle to DistanceHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DistanceHeight as text
%        str2double(get(hObject,'String')) returns contents of DistanceHeight as a double


% --- Executes during object creation, after setting all properties.
function DistanceHeight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DistanceHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupVec.
function popupVec_Callback(hObject, eventdata, handles)
% hObject    handle to popupVec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupVec contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupVec


% --- Executes during object creation, after setting all properties.
function popupVec_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupVec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in export.
function export_Callback(hObject, eventdata, handles)
E2 = handles.E2;
Dxf = handles.Dxf;
A = handles.A;
% dlmwrite('output.txt',E,'delimiter','\t', 'precision', '%.3f');
fileID = fopen('output.txt','w');
formatSpec = '%d %s %f %f %f %f %f \r\n';

[nrow, ncol] = size(E2);
for row = 1:nrow
    fprintf(fileID, formatSpec, E2{row,:});
end

D = Dxf(:,5:7);

Dx = D(:,1);
Dy = D(:,2);
Dz = D(:,3);

% Dx (row+1,1) = Dx (1,1);
% Dy (row+1,1) = Dy (1,1);
% Dz (row+1,1) = Dz (1,1);
% fileID = fopen('outputdxf.txt','w');
% formatSpec = '%f %f %f \r\n';
% 
% [nrow, ncol] = size(D);
% for row = 1:nrow
%     fprintf(fileID, formatSpec, D{row,:});
% end
Dx=cell2mat(Dx);
Dy=cell2mat(Dy);
Dz=cell2mat(Dz);
polydxf('outputdxf.dxf',Dx,Dy,Dz,A)
% hObject    handle to export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in clearDraw.
function clearDraw_Callback(hObject, eventdata, handles)
% if cam == 1
   axes(handles.Cam1);
   cla;
   index_selected = get(handles.listCam1,'Value');
    file_list = get(handles.listCam1,'String');
    filename = file_list{index_selected};
    % handles.index1 = index_selected;

    DirFile=handles.DirFile;
    complete = strcat(DirFile,'/',filename);

    J = imread(complete);
    axes(handles.Cam1);
    handles.Cam1 = imshow(J,[]);

    [~, file, ~] = fileparts (filename);
    file1 = file;
    save('file1.mat', 'file1');
% elseif cam == 2
   axes(handles.Cam2);
   cla;
   index_selected = get(handles.listCam2,'Value');
    file_list = get(handles.listCam2,'String');
    filename = file_list{index_selected};
    % handles.index1 = index_selected;

    DirFile=handles.DirFile;
    complete = strcat(DirFile,'/',filename);

    J = imread(complete);
    axes(handles.Cam2);
    handles.Cam2 = imshow(J,[]);

    [~, file, ~] = fileparts (filename);
    file2 = file;
    save('file2.mat', 'file2');
% X = str2num(get(handles.BuildingID,'String'));
% X = X+1;
% set(handles.BuildingID,'String',X);
% hObject    handle to clearDraw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in toggleVecDraw.
function toggleVecDraw_Callback(hObject, eventdata, handles)

i=1;
A=[];
while get(hObject,'Value')
   [x, y] = ginput(1);
   A(i,1) = x;
   A(i,2) = y;
   i = i+1;
end
B = A;
hold on
fill (B(1:i-3,1),B(1:i-3,2),'r')
hold on

storedStructure = load('Cam1.mat');
Cam1 = storedStructure.Cam1;

cam = get(handles.popupVec,'Value');

c = handles.c1;
x0 = handles.x01;
y0 = handles.y01;
K1 = handles.K1;
K2 = handles.K2;
K3 = handles.K3;
P1 = handles.P1;
P2 = handles.P2;

if cam == 1
    storedStructure2 = load('file1.mat');
    file1 = storedStructure2.file1;
    file1 = strcat(file1,'.JPG');
    row= Cam1(strcmp(file1,Cam1(:,1)),:);

    X0 = cell2mat(row (1,2));
    Y0 = cell2mat(row (1,3));
    Z0 = cell2mat(row (1,4));
    omega = cell2mat(row (1,5));
    phi = cell2mat(row (1,6));
    kappa = cell2mat(row (1,7));
    
    index_selected = get(handles.listCam1,'Value');
    file_list = get(handles.listCam1,'String');
    filename = file_list{index_selected};
elseif cam == 2
    storedStructure2 = load('file2.mat');
    file2 = storedStructure2.file2;
    file2 = strcat(file2,'.JPG');
    row= Cam1(strcmp(file2,Cam1(:,1)),:);

    X0 = cell2mat(row (1,2));
    Y0 = cell2mat(row (1,3));
    Z0 = cell2mat(row (1,4));
    omega = cell2mat(row (1,5));
    phi = cell2mat(row (1,6));
    kappa = cell2mat(row (1,7));
    
    index_selected = get(handles.listCam2,'Value');
    file_list = get(handles.listCam2,'String');
    filename = file_list{index_selected};
else
    disp('An error has occured')
end

pixw = handles.pixw;
pixh = handles.pixh;
% pixw = str2num(get(handles.pixwidth,'String'));
% pixh = str2num(get(handles.pixheight,'String'));
pix = str2num(get(handles.pixel,'String'));

ID = str2num(get(handles.BuildingID,'String'));

[ro,~] = size(B);
clic = ro - 2;

storedStructure2 = load('DSM.mat');
DSM = storedStructure2.DSM;

X0tiff = handles.X0tiff;
Y0tiff = handles.Y0tiff;
ResTiff = handles.ResTiff;


%bilinear interpolation for ground altitude 
%(1st initial value: X Y of camera)
% Z = 250;
Z = bilineinter (DSM,X0tiff,Y0tiff,ResTiff,X0,Y0);

r1=rotation2(omega,phi,kappa);

X0tiff = sprintf('%.4f',X0tiff);
X0tiff = str2double(X0tiff);

Y0tiff = sprintf('%.4f',Y0tiff);
Y0tiff = str2double(Y0tiff);

E={};
% realdot = @(u,v)u*transpose(v);
for i=1:ro-2
    
    xr = (B(i,1)*pix)-((pix*pixw)/2);
    yr = -(B(i,2)*pix)+((pix*pixh)/2);
    
    x = xr - x0;
    y = yr - y0;
    
    %calibration of image coordinates
    %radial corrections
    r = sqrt(x^2+y^2);
    % dr = K1*r^3-K2*r^5+K3*r^7;
    % drx = x1*(dr/r);
    % dry = y1*(dr/r);
    drx = x*(K1*r^2+K2*r^4+K3*r^6);
    dry = y*(K1*r^2+K2*r^4+K3*r^6);

    %decentering corrections
    ddx = P1*(3*x^2+y^2)+2*P2*x*y;
    ddy = P2*(x^2+3*y^2)+2*P1*x*y;

    x = xr -x0 - drx - ddx;
    y = yr -y0 - dry - ddy;

    for j=1:3
        %collinearity
        Xx = (Z-Z0)*((r1(1,1)*x + r1(2,1)*y - r1(3,1)*c)/(r1(1,3)*x + r1(2,3)*y - r1(3,3)*c));
        Xx = sprintf('%.4f',Xx);
        Xx = str2double(Xx);
        X0 = sprintf('%.4f',X0);
        X0 = str2double(X0);
        X = X0 + Xx;
    
        Yy = (Z-Z0)*((r1(1,2)*x + r1(2,2)*y - r1(3,2)*c)/(r1(1,3)*x + r1(2,3)*y - r1(3,3)*c));
        Yy = sprintf('%.4f',Yy);
        Yy = str2double(Yy);
        Y0 = sprintf('%.4f',Y0);
        Y0 = str2double(Y0);
        Y = Y0 + Yy;
    
        Z = bilineinter (DSM,X0tiff,Y0tiff,ResTiff,X,Y);
    end
    Xa = X;
    Ya = Y;
    Za = Z;
    
    xr2 = xr+((pix*pixw)/2);
    xr2 = xr2/pix;
    yr2 = yr+((pix*pixh)/2);
    yr2 = yr2/pix;
    
    E(i,1) = num2cell(i+(ID*100));
    E(i,2) = cellstr(filename);
    E(i,3) = num2cell(xr2);
    E(i,4) = num2cell(yr2);
    E(i,5) = num2cell(Xa);
    E(i,6) = num2cell(Ya);
    E(i,7) = num2cell(Za);
end

handles.E = E;
handles.clic = clic;
guidata(handles.toggleVecDraw, handles);




% hObject    handle to toggleVecDraw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of toggleVecDraw



function BuildingID_Callback(hObject, eventdata, handles)
% hObject    handle to BuildingID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BuildingID as text
%        str2double(get(hObject,'String')) returns contents of BuildingID as a double


% --- Executes during object creation, after setting all properties.
function BuildingID_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BuildingID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in BuildingIDAdd.
function BuildingIDAdd_Callback(hObject, eventdata, handles)
X = str2num(get(handles.BuildingID,'String'));
X = X+1;
set(handles.BuildingID,'String',X);
% hObject    handle to BuildingIDAdd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in BuildingIDSubs.
function BuildingIDSubs_Callback(hObject, eventdata, handles)
X = str2num(get(handles.BuildingID,'String'));
if X > 1
X = X-1;
else
    X=X;
end
set(handles.BuildingID,'String',X);
% hObject    handle to BuildingIDSubs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in buttonSaveDraw.
function buttonSaveDraw_Callback(hObject, eventdata, handles)
clic = handles.clic;
A = handles.A;
E = handles.E;
E2 = handles.E2;
Dxf = handles.Dxf;
E2 = [E2;E;{[] [] [] [] [] [] []}];

Dxf = [Dxf;E;E(1,1) E(1,2) E(1,3) E(1,4) E(1,5) E(1,6) E(1,7)];
A = [A;clic];

handles.A = A;
handles.E2 = E2;
handles.Dxf = Dxf;
guidata(handles.buttonSaveDraw, handles);
% hObject    handle to buttonSaveDraw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in PointSelect2.
function PointSelect2_Callback(hObject, eventdata, handles)
set(handles.zoomto,'Enable', 'On');
[x1,y1] = ginput (1);
handles.X = x1;
handles.Y = y1;
hold on;
plot(x1,y1,'y+');
hold on;
guidata(handles.PointSelect,handles);
pixw = handles.pixw;
pixh = handles.pixh;
% pixw = str2num(get(handles.pixwidth,'String'));
% pixh = str2num(get(handles.pixheight,'String'));
pix = str2num(get(handles.pixel,'String'));

x1=x1*pix;
x1=x1-((pix*pixw)/2);
y1=y1*pix;
y1=-y1+((pix*pixh)/2);

storedStructure = load('Cam1.mat');
Cam1 = storedStructure.Cam1;

% hObject    handle to PointSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%interior orientation
c1 = handles.c1;
x01 = handles.x01;
y01 = handles.y01;
K1 = handles.K1;
K2 = handles.K2;
K3 = handles.K3;
P1 = handles.P1;
P2 = handles.P2;


%exterior orientation
storedStructure2 = load('file2.mat');
file2 = storedStructure2.file2;
file2 = strcat(file2,'.JPG');
 row = Cam1(strcmp(file2,Cam1(:,1)),:);
% row = find (ismember(Cam1,file1),1)

X01 = cell2mat(row (1,2));
Y01 = cell2mat(row (1,3));
Z01 = cell2mat(row (1,4));
omega1 = cell2mat(row (1,5));
phi1 = cell2mat(row (1,6));
kappa1 = cell2mat(row (1,7));

%rotation matrix
r1=rotation2(omega1,phi1,kappa1);

%bilinear interpolation for ground altitude 
%(1st initial value: X Y of camera)
% Z = 250;
storedStructure2 = load('DSM.mat');
DSM = storedStructure2.DSM;

X0tiff = handles.X0tiff;
Y0tiff = handles.Y0tiff;
ResTiff = handles.ResTiff;

Z = bilineinter (DSM,X0tiff,Y0tiff,ResTiff,X01,Y01);

x = x1-x01;
y = y1-y01;

%calibration of image coordinates
%radial corrections
r = sqrt(x^2+y^2);
% dr = K1*r^3-K2*r^5+K3*r^7;
% drx = x1*(dr/r);
% dry = y1*(dr/r);
drx = x*(K1*r^2+K2*r^4+K3*r^6);
dry = y*(K1*r^2+K2*r^4+K3*r^6);

%decentering corrections
ddx = P1*(3*x^2+y^2)+2*P2*x*y;
ddy = P2*(x^2+3*y^2)+2*P1*x*y;

x = x - drx - ddx;
y = y - dry - ddy;

%loop to iterate Z
for i=1:3
    %collinearity
    Xx = (Z-Z01)*((r1(1,1)*x + r1(2,1)*y - r1(3,1)*c1)/(r1(1,3)*x + r1(2,3)*y - r1(3,3)*c1));
    Xx = sprintf('%.4f',Xx);
    Xx = str2double(Xx);
    X01 = sprintf('%.4f',X01);
    X01 = str2double(X01);
    X = X01 + Xx;
    
    Yy = (Z-Z01)*((r1(1,2)*x + r1(2,2)*y - r1(3,2)*c1)/(r1(1,3)*x + r1(2,3)*y - r1(3,3)*c1));
    Yy = sprintf('%.4f',Yy);
    Yy = str2double(Yy);
    Y01 = sprintf('%.4f',Y01);
    Y01 = str2double(Y01);
    Y = Y01 + Yy;
    
    Z = bilineinter (DSM,X0tiff,Y0tiff,ResTiff,X,Y);
end

set (handles.editXres, 'string', sprintf('%0.3f',X));
set (handles.editYres, 'string', sprintf('%0.3f',Y));
set (handles.editZres, 'string', sprintf('%0.3f',Z));

storedStructure2 = load('file1.mat');
    file1 = storedStructure2.file1;
    file1 = strcat(file1,'.JPG');
    row= Cam1(strcmp(file1,Cam1(:,1)),:);

    X02 = cell2mat(row (1,2));
    Y02 = cell2mat(row (1,3));
    Z02 = cell2mat(row (1,4));
    omega2 = cell2mat(row (1,5));
    phi2 = cell2mat(row (1,6));
    kappa2 = cell2mat(row (1,7));
    
    
c2 = handles.c1;
x02 = handles.x01;
y02 = handles.y01;

%rotation matrix
r2 = rotation2(omega2, phi2, kappa2);

%collinearity
Xt2 = X-X02;
Yt2 = Y-Y02;
Zt2 = Z-Z02;

x2 = x02 - c2*((r2(1,1)*Xt2 + r2(1,2)*Yt2 + r2(1,3)*Zt2)/(r2(3,1)*Xt2 + r2(3,2)*Yt2 + r2(3,3)*Zt2)) ;
y2 = y02 - c2*((r2(2,1)*Xt2 + r2(2,2)*Yt2 + r2(2,3)*Zt2)/(r2(3,1)*Xt2 + r2(3,2)*Yt2 + r2(3,3)*Zt2)) ;

%calibration of image coordinates
%radial corrections
r = sqrt(x2^2+y2^2);
% dr = K1*r^3-K2*r^5+K3*r^7;
% drx = x1*(dr/r);
% dry = y1*(dr/r);
drx = x2*(K1*r^2+K2*r^4+K3*r^6);
dry = y2*(K1*r^2+K2*r^4+K3*r^6);

%decentering corrections
ddx = P1*(3*x2^2+y2^2)+2*P2*x2*y2;
ddy = P2*(x2^2+3*y2^2)+2*P1*x2*y2;

x2 = x2 + drx + ddx;
y2 = y2 + dry + ddy;


x2 = x2+((pix*pixw)/2);
x2 = x2/pix;
y2 = -(y2-((pix*pixh)/2));
y2 = y2/pix;

axes(handles.Cam1)
hold on;
plot(x2,y2,'rx');
hold on;
% hObject    handle to PointSelect2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
