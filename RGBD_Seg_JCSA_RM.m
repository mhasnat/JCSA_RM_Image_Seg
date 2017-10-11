function varargout = RGBD_Seg_JCSA_RM(varargin)
% RGBD_SEG_JCSA_RM MATLAB code for RGBD_Sef_JCSA_RM.fig
%      RGBD_SEG_JCSA_RM, by itself, creates a new RGBD_SEG_JCSA_RM or raises the existing
%      singleton*.
%
%      H = RGBD_SEG_JCSA_RM returns the handle to a new RGBD_SEG_JCSA_RM or the handle to
%      the existing singleton*.
%
%      RGBD_SEG_JCSA_RM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RGBD_SEG_JCSA_RM.M with the given input arguments.
%
%      RGBD_SEG_JCSA_RM('Property','Value',...) creates a new RGBD_SEG_JCSA_RM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RGBD_Seg_JCSA_RM_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RGBD_Seg_JCSA_RM_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RGBD_Seg_JCSA_RM

% Last Modified by GUIDE v2.5 11-Oct-2017 11:31:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @RGBD_Seg_JCSA_RM_OpeningFcn, ...
    'gui_OutputFcn',  @RGBD_Seg_JCSA_RM_OutputFcn, ...
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


% --- Executes just before RGBD_Seg_JCSA_RM is made visible.
function RGBD_Seg_JCSA_RM_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RGBD_Seg_JCSA_RM (see VARARGIN)

% Choose default command line output for RGBD_Seg_JCSA_RM
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes RGBD_Seg_JCSA_RM wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = RGBD_Seg_JCSA_RM_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in save_seg.
function save_seg_Callback(hObject, eventdata, handles)
% hObject    handle to save_seg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Select file name
FileName = uiputfile('*.*');

if(~FileName)
    return;
end

imwrite(handles.segres, strcat(FileName, '.png'));

% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
warning off
addpath('exportfig');
addpath('vmfmatlab');
addpath('rgbd');
handles.MethodType = 'rb_jcsa_rm';

guidata(hObject,handles);

% --- Executes on button press in cluster_normals.
function cluster_normals_Callback(hObject, eventdata, handles)
% hObject    handle to cluster_normals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
display('Segmenting Image ... please wait ...')
% Get the radiobutton informations
if(get(handles.rb_sc_full,'Value'))
    sc = 1;
elseif(get(handles.rb_sc_quarter,'Value'))
    sc = 4;
else
    sc = 2;
end

% Load the threshold values
kMax = str2num(get(handles.tb_k_max, 'String'));

thOptions.thDivNormalMax = str2num(get(handles.th_dist_region, 'String'));
thOptions.thDivNormalMin = 1;
thOptions.planarityTh = str2num(get(handles.th_planar, 'String'));
thOptions.thKappa = str2num(get(handles.th_kappa_pl, 'String'));
thOptions.edgeStrengthTh = str2num(get(handles.th_boundary, 'String'));;

opt.showLLH = 0;
opt.showIt = 0;
opt.numiter = 20;

% Load necessary data from global data
handles = guidata(hObject);
allInfo = handles.allInfo;

szh(1) = uint16(size(handles.rgbImg(1:sc:end,1:sc:end,:),1));
szh(2) = uint16(size(handles.rgbImg(1:sc:end,1:sc:end,:),2));

% accumulate features
combFeat = getCombinedFeatures(allInfo, sc);

% Generate oversegmentation
if(strcmp(handles.MethodType, 'rb_jcsa_rm') | strcmp(handles.MethodType, 'rb_jcsa'))
    display('Applying JCSA Clustering ...');
    label = uint8(fusionBD_Color_3D_Axis(combFeat, kMax, [1 1 1], opt));
elseif(strcmp(handles.MethodType, 'rb_jcsd_rm') | strcmp(handles.MethodType, 'rb_jcsd'))
    display('Applying JCSD Clustering ...');
    label = uint8(fusionBD_Color_3D_Normal(combFeat, kMax, [1 1 1], opt));
else
    display('Choose your method name properly');
end

labelImg = reshape(label, szh(1), szh(2));

% post-processing
display('Applying Post-processing ...');
bdry = seg2bdry_2(labelImg);
regImg = bwlabeln(~bdry);
rgbImg = allInfo.rgbImg(1:2:end, 1:2:end, :);

[edges, tmp, neighbors, tmpLabImg] = seg2fragments(double(regImg), rgbImg, 10);

display('Applying Region-Megring method ...');
if(strcmp(handles.MethodType, 'rb_jcsa_rm'))
    % Compute parameters of the distributions
    parNormal = getKappaWMM(combFeat(:,7:9), tmpLabImg(:));
    parColor = getParamsDivergenceGMM(combFeat(:,1:3), tmpLabImg(:));
    
    % Set initial parameters for the neighborhood regions
    rgnb.nbsegs = getRegionNeighbors(tmpLabImg);
    rgnb.Div_N = ones(length(rgnb.nbsegs))*5000;
    rgnb.Div_C = ones(length(rgnb.nbsegs))*5000;
    rgnb.KappaMerged = ones(length(rgnb.nbsegs)) * -20;
    
    % Update parameters for the neighborhood regions
    for i=1:length(rgnb.nbsegs)
        adjs = rgnb.nbsegs{i};
        adjEdges = cell(length(adjs),1);
        
        for j = 1:length(adjs)
            commonEdge = intersect(neighbors.segment_fragments{i}, neighbors.segment_fragments{adjs(j)});
            commIndx = [];
            for k=1:length(commonEdge)
                commIndx = [commIndx; fliplr(floor(edges{commonEdge(k)}))];
            end
            rgnb.adjEdges{i,adjs(j)} = commIndx;
            
            % Compute kappa after it is merged with the neighbor cluster
            rgnb.KappaMerged(i,adjs(j)) = mergeKappaWMM(parNormal, [i,adjs(j)]);
            
            % Compute BD_normal among the neighbors
            rgnb.Div_N(i,adjs(j)) = parNormal.DLNF(i) - parNormal.DLNF(adjs(j)) - ( (parNormal.eta(i,:) - parNormal.eta(adjs(j),:)) * parNormal.theta_cl(adjs(j),:)');
            
            % Compute BD_color among the neighbors
            rgnb.Div_C(i,adjs(j)) = compute_Div_C(parColor, i, adjs(j));
        end
    end
    
    % Apply RM method
    img = getRAGMergeSegments(combFeat, tmpLabImg, thOptions, parNormal, parColor, rgnb, allInfo, sc, handles.MethodType);
elseif(strcmp(handles.MethodType, 'rb_jcsd_rm'))
    % Compute parameters of the distributions
    parNormal = getParamsVMFBD(combFeat(:,7:9), tmpLabImg(:));
    parColor = getParamsDivergenceGMM(combFeat(:,1:3), tmpLabImg(:));
    
    % Set initial parameters for the neighborhood regions
    rgnb.nbsegs = getRegionNeighbors(tmpLabImg);
    rgnb.Div_N = ones(length(rgnb.nbsegs))*5000;
    rgnb.Div_C = ones(length(rgnb.nbsegs))*5000;
    rgnb.KappaMerged = ones(length(rgnb.nbsegs)) * -20;
    
    % Update parameters for the neighborhood regions
    for i=1:length(rgnb.nbsegs)
        adjs = rgnb.nbsegs{i};
        adjEdges = cell(length(adjs),1);
        
        for j = 1:length(adjs)
            commonEdge = intersect(neighbors.segment_fragments{i}, neighbors.segment_fragments{adjs(j)});
            commIndx = [];
            for k=1:length(commonEdge)
                commIndx = [commIndx; fliplr(floor(edges{commonEdge(k)}))];
            end
            rgnb.adjEdges{i,adjs(j)} = commIndx;
            
            % Compute kappa after it is merged with the neighbor cluster
            rgnb.KappaMerged(i,adjs(j)) = mergeKappaVMFMM(parNormal, [i,adjs(j)]);
            
            % Compute BD_normal among the neighbors
            rgnb.Div_N(i,adjs(j)) = parNormal.DLNF(i) - parNormal.DLNF(adjs(j)) - ( (parNormal.eta(i,:) - parNormal.eta(adjs(j),:)) * parNormal.theta_cl(adjs(j),:)');
            
            % Compute BD_color among the neighbors
            rgnb.Div_C(i,adjs(j)) = compute_Div_C(parColor, i, adjs(j));
        end
    end
    
    % Apply RM method
    img = getRAGMergeSegments(combFeat, tmpLabImg, thOptions, parNormal, parColor, rgnb, allInfo, sc, handles.MethodType);
else
    img = tmpLabImg;
end

display('Displaying segmentation results ...');
handles.segres = label2rgb(assignRandomLabel(img));
guidata(hObject,handles);

imshow(handles.segres, 'parent',handles.ax_segres);

%%

% % Disable clustering button
% set(handles.cluster_normals,'Enable','off');


% --- Executes on button press in rb_allInfoImg.
function rb_allInfoImg_Callback(hObject, eventdata, handles)
% hObject    handle to rb_allInfoImg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb_allInfoImg

% Set radio buttons for normals as on
set(handles.rb_allInfoImg,'Value',1);

% Set all other radio buttons off
set(handles.rb_colimg,'Value',0);
set(handles.rb_img_normal,'Value',0);
set(handles.rb_depImg,'Value',0);

% --- Executes on button press in rb_img_normal.
function rb_img_normal_Callback(hObject, eventdata, handles)
% hObject    handle to rb_img_normal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb_img_normal

% Set radio buttons for normals as on
set(handles.rb_img_normal,'Value',1);

% Set all other radio buttons off
set(handles.rb_allInfoImg,'Value',0);
set(handles.rb_depImg,'Value',0);
set(handles.rb_colimg,'Value',0);

% --- Executes on button press in rb_depImg.
function rb_depImg_Callback(hObject, eventdata, handles)
% hObject    handle to rb_depImg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb_depImg

% Set radio buttons for depth as on
set(handles.rb_depImg,'Value',1);

% Set all other radio buttons off
set(handles.rb_allInfoImg,'Value',0);
set(handles.rb_img_normal,'Value',0);
set(handles.rb_colimg,'Value',0);

% --- Executes on button press in rb_colImg.
function rb_colImg_Callback(hObject, eventdata, handles)
% hObject    handle to rb_colImg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb_colImg

% Set radio buttons for color as on
set(handles.rb_colimg,'Value',1);

% Set all other radio buttons off
set(handles.rb_img_normal,'Value',0);
set(handles.rb_allInfoImg,'Value',0);
set(handles.rb_depImg,'Value',0);


% --- Executes on button press in rb_colimg.
function rb_colimg_Callback(hObject, eventdata, handles)
% hObject    handle to rb_colimg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb_colimg

% Set radio buttons for normals as on
set(handles.rb_colimg,'Value',1);

% Set all other radio buttons off
set(handles.rb_allInfoImg,'Value',0);
set(handles.rb_img_normal,'Value',0);
set(handles.rb_depImg,'Value',0);


% --- Executes on button press in pb_load_image.
function pb_load_image_Callback(hObject, eventdata, handles)
% hObject    handle to pb_load_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Select the appropriate file
[FileName,PathName,FilterIndex] = uigetfile;

if(~FileName)
    return;
end

% Load samples and labels
load(strcat(PathName, FileName));

handles.rgbImg = rgbd_data.rgbImg;
handles.depImg = rgbd_data.depImg;
handles.imgNormals = rgbd_data.imgNormals;

handles.allInfo = rgbd_data;

guidata(hObject,handles);

% Display images
imshow(handles.rgbImg, 'parent',handles.ax_colImg);
imshow(handles.depImg, [],'parent',handles.ax_depthImg);
imshow(handles.imgNormals, 'parent',handles.ax_normal);

% Enable clustering button
set(handles.cluster_normals,'Enable','on');


% --- Executes on button press in rb_sc_half.
function rb_sc_half_Callback(hObject, eventdata, handles)
% hObject    handle to rb_sc_half (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb_sc_half

% Set radio buttons for normals as on
set(handles.rb_sc_half,'Value',1);

% Set all other radio buttons off
set(handles.rb_sc_full,'Value',0);
set(handles.rb_sc_quarter,'Value',0);


% --- Executes on button press in rb_sc_full.
function rb_sc_full_Callback(hObject, eventdata, handles)
% hObject    handle to rb_sc_full (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb_sc_full

% Set radio buttons for normals as on
set(handles.rb_sc_full,'Value',1);

% Set all other radio buttons off
set(handles.rb_sc_half,'Value',0);
set(handles.rb_sc_quarter,'Value',0);

% --- Executes on button press in rb_sc_quarter.
function rb_sc_quarter_Callback(hObject, eventdata, handles)
% hObject    handle to rb_sc_quarter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb_sc_quarter

% Set radio buttons for normals as on
set(handles.rb_sc_quarter,'Value',1);

% Set all other radio buttons off
set(handles.rb_sc_half,'Value',0);
set(handles.rb_sc_full,'Value',0);


% --- Executes on button press in rb_nc_auto.
function rb_nc_auto_Callback(hObject, eventdata, handles)
% hObject    handle to rb_nc_auto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb_nc_auto

% Set radio buttons for normals as on
set(handles.rb_nc_auto,'Value',1);

% Set all other radio buttons off
set(handles.rb_nc_manual,'Value',0);
set(handles.rbg_auto_num_comp, 'Visible', 'on');
set(handles.ax_res_analysis, 'Visible', 'on');
set(handles.tb_numcomp, 'Visible', 'off');

handles.compselectiontype = 1; % means automatic component selection
guidata(hObject,handles);

% --- Executes on button press in rb_nc_manual.
function rb_nc_manual_Callback(hObject, eventdata, handles)
% hObject    handle to rb_nc_manual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb_nc_manual

% Set radio buttons for normals as on
set(handles.rb_nc_manual,'Value',1);

% Set all other radio buttons off
set(handles.rb_nc_auto,'Value',0);
set(handles.tb_numcomp, 'Visible', 'on');
set(handles.rbg_auto_num_comp, 'Visible', 'off');
cla(handles.ax_res_analysis,'reset');
set(handles.ax_res_analysis, 'Visible', 'off');

handles.compselectiontype = 0; % means number of component set externally
guidata(hObject,handles);

% --- Executes when selected object is changed in rbg_auto_num_comp.
function rbg_auto_num_comp_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in rbg_auto_num_comp
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

h_selectedRadioButton = get(handles.rbg_auto_num_comp,'SelectedObject');
handles.selectedRadioTag = get(h_selectedRadioButton,'tag');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function rbg_auto_num_comp_CreateFcn(hObject, eventdata, handles)

function tb_numcomp_Callback(hObject, eventdata, handles)

function tb_numcomp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in rbg_seg_method.
function rbg_seg_method_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in rbg_seg_method
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

h_selectedRadioButton = get(handles.rbg_seg_method,'SelectedObject');
handles.MethodType = get(h_selectedRadioButton,'tag');
guidata(hObject,handles);



function tb_k_max_Callback(hObject, eventdata, handles)
% hObject    handle to tb_k_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tb_k_max as text
%        str2double(get(hObject,'String')) returns contents of tb_k_max as a double


% --- Executes during object creation, after setting all properties.
function tb_k_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tb_k_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function th_kappa_pl_Callback(hObject, eventdata, handles)
% hObject    handle to th_kappa_pl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of th_kappa_pl as text
%        str2double(get(hObject,'String')) returns contents of th_kappa_pl as a double


% --- Executes during object creation, after setting all properties.
function th_kappa_pl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to th_kappa_pl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function th_dist_region_Callback(hObject, eventdata, handles)
% hObject    handle to th_dist_region (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of th_dist_region as text
%        str2double(get(hObject,'String')) returns contents of th_dist_region as a double


% --- Executes during object creation, after setting all properties.
function th_dist_region_CreateFcn(hObject, eventdata, handles)
% hObject    handle to th_dist_region (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function th_boundary_Callback(hObject, eventdata, handles)
% hObject    handle to th_boundary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of th_boundary as text
%        str2double(get(hObject,'String')) returns contents of th_boundary as a double


% --- Executes during object creation, after setting all properties.
function th_boundary_CreateFcn(hObject, eventdata, handles)
% hObject    handle to th_boundary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function th_planar_Callback(hObject, eventdata, handles)
% hObject    handle to th_planar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of th_planar as text
%        str2double(get(hObject,'String')) returns contents of th_planar as a double


% --- Executes during object creation, after setting all properties.
function th_planar_CreateFcn(hObject, eventdata, handles)
% hObject    handle to th_planar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in rb_jcsa.
function rb_jcsa_Callback(hObject, eventdata, handles)
% hObject    handle to rb_jcsa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb_jcsa
handles.MethodType = 'rb_jcsa';


% --- Executes on button press in rb_jcsd.
function rb_jcsd_Callback(hObject, eventdata, handles)
% hObject    handle to rb_jcsd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb_jcsd
handles.MethodType = 'rb_jcsd';


% --- Executes on button press in rb_jcsa_rm.
function rb_jcsa_rm_Callback(hObject, eventdata, handles)
% hObject    handle to rb_jcsa_rm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb_jcsa_rm
handles.MethodType = 'rb_jcsa_rm';


% --- Executes on button press in rb_jcsd_rm.
function rb_jcsd_rm_Callback(hObject, eventdata, handles)
% hObject    handle to rb_jcsd_rm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb_jcsd_rm
handles.MethodType = 'rb_jcsd_rm';
