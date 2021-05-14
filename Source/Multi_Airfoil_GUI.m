function varargout = Multi_Airfoil_GUI(varargin)
% MULTI_AIRFOIL_GUI MATLAB code for Multi_Airfoil_GUI.fig
%      MULTI_AIRFOIL_GUI, by itself, creates a new MULTI_AIRFOIL_GUI or raises the existing
%      singleton*.
%
%      H = MULTI_AIRFOIL_GUI returns the handle to a new MULTI_AIRFOIL_GUI or the handle to
%      the existing singleton*.
%
%      MULTI_AIRFOIL_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MULTI_AIRFOIL_GUI.M with the given input arguments.
%
%      MULTI_AIRFOIL_GUI('Property','Value',...) creates a new MULTI_AIRFOIL_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Multi_Airfoil_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Multi_Airfoil_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Multi_Airfoil_GUI

% Last Modified by GUIDE v2.5 27-Apr-2021 14:33:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Multi_Airfoil_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @Multi_Airfoil_GUI_OutputFcn, ...
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


% --- Executes just before Multi_Airfoil_GUI is made visible.
function Multi_Airfoil_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Multi_Airfoil_GUI (see VARARGIN)

% Choose default command line output for Multi_Airfoil_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%Initialize
evalin('base','clear all');
evalin('base','clc');

% UIWAIT makes Multi_Airfoil_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Multi_Airfoil_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




function input1_Callback(hObject, eventdata, handles)
% hObject    handle to input1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input1 as text
%        str2double(get(hObject,'String')) returns contents of input1 as a double


% --- Executes during object creation, after setting all properties.
function input1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function Dis1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Dis1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Dis2_Callback(hObject, eventdata, handles)
% hObject    handle to Dis2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Dis2 as text
%        str2double(get(hObject,'String')) returns contents of Dis2 as a double


% --- Executes during object creation, after setting all properties.
function Dis2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Dis2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function input2_Callback(hObject, eventdata, handles)
% hObject    handle to input2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input2 as text
%        str2double(get(hObject,'String')) returns contents of input2 as a double


% --- Executes during object creation, after setting all properties.
function input2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function input3_Callback(hObject, eventdata, handles)
% hObject    handle to input3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input3 as text
%        str2double(get(hObject,'String')) returns contents of input3 as a double


% --- Executes during object creation, after setting all properties.
function input3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function input4_Callback(hObject, eventdata, handles)
% hObject    handle to input4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input4 as text
%        str2double(get(hObject,'String')) returns contents of input4 as a double


% --- Executes during object creation, after setting all properties.
function input4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Dis1_Callback(hObject, eventdata, handles)
% hObject    handle to Dis1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Dis1 as text
%        str2double(get(hObject,'String')) returns contents of Dis1 as a double


% --- Executes on button press in download_button1.
function download_button1_Callback(hObject, eventdata, handles)
% hObject    handle to download_button1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Plot_type title load_success plot_success
global AF_name numPtsAF numAF 
global XX YY XC XB_AF YB_AF
global Vx Vy Cp CpXY 
Plot_type = 'Airfoil_Dat';
[plot_success] = output_gui(Plot_type,XB_AF,YB_AF,XC,Cp,XX,YY,Vx,Vy,CpXY,AF_name,numPtsAF,numAF,title);

if plot_success == 1 & load_success == 1
    questdlg('Airfoil Dat. Download Complete',...
            'Confirm',...
            'OK','OK');
else
    questdlg('Download Fail',...
            'Confirm',...
            'OK','OK');
end


% --- Executes on button press in download_button2.
function download_button2_Callback(hObject, eventdata, handles)
% hObject    handle to download_button2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Plot_type title load_success cal_success
global AF_name numPtsAF numAF 
global XX YY XC XB_AF YB_AF
global Vx Vy Cp CpXY  
Str = get(handles.popupmenu1 , 'String');
Val = get(handles.popupmenu1 , 'Value');
Plot_type = Str{Val}; 
[cal_success] = output_gui(Plot_type,XB_AF,YB_AF,XC,Cp,XX,YY,Vx,Vy,CpXY,AF_name,numPtsAF,numAF,title);

if cal_success == 1 & load_success == 1
    questdlg('Plot Dat. Download Complete',...
            'Confirm',...
            'OK','OK');
else
    questdlg('Download Fail',...
            'Confirm',...
            'OK','OK');
end

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Set Global Variable
global Plot_type title load_success 
global AF_name numAF numPtsAF panInd 
global xVals yVals stepsize maxVert 
global Ysl
global XC XB YB XX YY XYB XB_AF YB_AF 
global Vx Vy Cp CpXY

set(handles.Dis1,'String',sprintf('\n\n\nLoading...'));
cla(handles.Dis2,'reset');
cla(handles.Dis3,'reset');

% Condition
AF_name = {get(handles.input1 , 'String')};
AF_load   = [1];
AF_flip   = [1];                                                            % Flip the airfoil vertically [1 = NO, -1 = YES]
AF_scale  = [1];                                                            % Scale airfoil
AF_angle  = [str2double(get(handles.input3 , 'String'))];                   % Angle of rotation of airfoil [deg]
AF_offset = [0 0];                                                          % Offset from origin [x y]

% User-defined knowns
Vinf = str2double(get(handles.input2 , 'String'));                          % Freestream velocity []  (just leave this at 1)

%Check Input Validation
set(handles.Dis1,'String','Loading....'); 
index = size(get(handles.Dis1,'string'), 1);                            % Get how many items are in the list box
set(handles.Dis1,'Value',index);
drawnow; 
[XB,YB,numPtsAF,numPanAF,title,string,load_success] = check_input_gui(AF_name,AF_load,AF_flip,AF_scale,AF_angle,AF_offset,Vinf);
    
if load_success == 0
    set(handles.Dis1,'String',string); 
    index = size(get(handles.Dis1,'string'), 1);                            % Get how many items are in the list box
    set(handles.Dis1,'Value',index);
    drawnow;
    return;
elseif load_success == -1
    set(handles.Dis1,'String',string); 
    index = size(get(handles.Dis1,'string'), 1);                            % Get how many items are in the list box
    set(handles.Dis1,'Value',index);
    drawnow;
    return;
end

set(handles.Dis1,'String',string); 
index = size(get(handles.Dis1,'String'), 1);                                % Get how many items are in the list box
set(handles.Dis1,'Value',index);
drawnow;
    %% Main Calculation 
    nGridX = str2double(get(handles.input4 , 'String'));
    nGridY = nGridX;

    % KNOWNS
    AoA  = 0;                                                               % Default Angle of attack [deg]
    
    % Number of airfoils
    numAF = length(AF_name);                                                % Number of airfoils
    
    for i=1:numAF
        title = char(strcat(title,AF_name(i),'_'));
    end
    
    %=========================================================================%
    %            Get Airfoil GEO. from check_input_gui.m                      %
    %=========================================================================%

    %% PANEL CALCULATIONS
    totPan   = sum(numPtsAF)-1;                                                 % Total number of panels (including false panels)
    falsePan = cumsum(numPanAF) + (1:1:numAF)';                                 % Indices of false panels
    falsePan = falsePan(1:end-1);                                               % Get rid of last entry in array

    % Array of actual airfoil panel indices
    panInd           = (1:1:totPan)';                                           % List of all panels
    panInd(falsePan) = 0;                                                       % Set intermediate panels to a value of zero
    panInd           = panInd(panInd ~= 0);                                     % Panel indices of the actual panels (not intermediate panels)
    iInd             = panInd;                                                  % Rename for use in functions
    jInd             = panInd;                                                  % Rename for use in functions

    % Starting/ending indices for referencing each airfoil's indices in panInd
    kS(1,1) = 1;                                                                % Set starting index for first airfoil
    for k = 1:1:numAF                                                           % Loop over all airfoils
        if (k == 1)                                                             % If it's the first airfoil
            kE(k,1) = kS(k) + numPanAF(k) - 1;                                  % Ending index of the first airfoil
        else                                                                    % If it's not the first airfoil
            kS(k,1) = kS(k-1) + numPanAF(k-1);                                  % Starting index of the rest of the airfoils
            kE(k,1) = kS(k) + numPanAF(k) - 1;                                  % Ending index of the rest of the airfoils
        end
    end

    %% CHECK PANEL DIRECTIONS - FLIP IF NECESSARY

    % Check for direction of points on each airfoil (We want Panels are CW) 
    % - Flips only airfoils that are incorrectly oriented
    for i = 1:1:numAF                                                           % Loop over all airfoils
        edge = zeros(numPanAF(i),1);                                            % Initialize edge value array
        for j = 1:1:numPanAF(i)                                                 % Loop over all panels
            edge(j,1) = (XB{i}(j+1)-XB{i}(j))*(YB{i}(j+1)+YB{i}(j));            % Compute edge values
        end
        sumEdge = sum(edge);                                                    % Sum all edge values

        if (sumEdge < 0)                                                        % If panels are CCW
            XB{i} = flipud(XB{i});                                              % Flip the X-data array
            YB{i} = flipud(YB{i});                                              % Flip the Y-data array
        end
    end

    % Individual airfoils to save for later plotting
    XB_AF = XB;                                                                 % X boundary points for each airfoil
    YB_AF = YB;                                                                 % Y boundary points for each airfoil

    % Consolidate all airfoils into one array
    XBT = [];                                                                   % Initialize temporary array
    YBT = [];                                                                   % Initialize temporary array
    for i = 1:1:numAF                                                           % Loop over all airfoils
        XBT = [XBT; XB{i}];                                                     % Concatenate all airfoil X boundary points
        YBT = [YBT; YB{i}];                                                     % Concatenate all airfoil Y boundary points
    end
    XB = XBT;                                                                   % Overwrite X boundary points
    YB = YBT;                                                                   % Overwrite Y boundary points

    % Number of panels
    numPtsTot = length(XB);                                                     % Number of boundary/control points
    numPanTot = length(XB)-1;                                                   % Number of panels

    %% PANEL METHOD GEOMETRY - REF [1]

    % Initialize variables
    XC    = zeros(numPanTot,1);                                                 % Initialize control point X-coordinate array
    YC    = zeros(numPanTot,1);                                                 % Initialize control point Y-coordinate array
    S     = zeros(numPanTot,1);                                                 % Intialize panel length array
    phiD  = zeros(numPanTot,1);                                                 % Initialize panel orientation angle array [deg]

    % Find geometric quantities of the airfoils
    for i = 1:1:numPanTot                                                       % Loop over all panels
        XC(i)   = 0.5*(XB(i)+XB(i+1));                                          % X-value of control point
        YC(i)   = 0.5*(YB(i)+YB(i+1));                                          % Y-value of control point
        dx      = XB(i+1)-XB(i);                                                % Change in X between boundary points
        dy      = YB(i+1)-YB(i);                                                % Change in Y between boundary points
        S(i)    = sqrt(dx^2 + dy^2);                                            % Length of the panel
        phiD(i) = atan2d(dy,dx);                                                % Angle of the panel (positive X-axis to inside face) [deg]
        if (phiD(i) < 0)                                                        % If panel angles are negative [deg]
            phiD(i) = phiD(i) + 360;                                            % Make all panel angles positive [deg]
        end
    end

    % Compute angle of panel normal w.r.t horizontal and include AoA
    deltaD             = phiD + 90;                                             % Angle from positive X-axis to outward normal vector [deg]
    betaD              = deltaD - AoA;                                          % Angle between freestream vector and outward normal vector [deg]
    betaD(betaD > 360) = betaD(betaD > 360) - 360;                              % Make all panel angles between 0 and 360 [deg]

    % Convert angles from [deg] to [rad]
    phi  = phiD.*(pi/180);                                                      % Convert from [deg] to [rad]
    beta = betaD.*(pi/180);                                                     % Convert from [deg] to [rad]

    % Exclude panels that are not on airfoils
    beta_Real = beta(panInd);                                                   % Keep values that are part of actual airfoils
    S_Real    = S(panInd);                                                      % Keep values that are part of actual airfoils

    %% COMPUTE SOURCE AND VORTEX PANEL STRENGTHS - REFS [2,3,6,7,10]

    % Number of panels
    numPan = sum(numPanAF);                                                     % Number of real airfoil panels

    % Geometric integral
    [I,J] = COMPUTE_IJ_SPM_N(XC,YC,XB,YB,phi,S,numPan,iInd,jInd);               % Call COMPUTE_IJ_SPM function (Refs [2] and [3])
    [K,L] = COMPUTE_KL_VPM_N(XC,YC,XB,YB,phi,S,numPan,iInd,jInd);               % Call COMPUTE_KL_VPM function (Refs [6] and [7])

    % Populate A matrix
    % - Simpler option: A = I + pi*eye(numPan,numPan);
    A = zeros(numPan,numPan);                                                   % Initialize the A matrix
    for i = 1:1:numPan                                                          % Loop over all i panels
        for j = 1:1:numPan                                                      % Loop over all j panels
            if (i == j)                                                         % If the panels are the same
                A(i,j) = pi;                                                    % Set A equal to pi
            else                                                                % If panels are not the same
                A(i,j) = I(i,j);                                                % Set A equal to I
            end
        end
    end

    % Right column of A matrix
    for k = 1:1:numAF                                                           % Loop over all airfoils
        for i = 1:1:numPan                                                      % Loop over all i panels (rows)
            A(i,numPan+k) = -sum(K(i,kS(k):kE(k)));                             % Add gamma term to right-most column of A matrix
        end
    end

    % Populate b array
    b = zeros(numPan,1);                                                        % Initialize the b array
    for i = 1:1:numPan                                                          % Loop over all i panels (rows)
        b(i) = -Vinf*2*pi*cos(beta_Real(i));                                    % Compute RHS array
    end

    % ===== Start/Enforce the Kutta Condition, REF [10] =====
    % \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

    % Bottom row of A matrix (Kutta condition)
    for k = 1:1:numAF                                                           % Loop over all airfoils
        for j = 1:1:numPan                                                      % Loop over all j panels (columns)
            A(numPan+k,j) = J(kS(k),j) + J(kE(k),j);                            % Source panel terms for each Kutta condition equation in A matrix
        end
        A(numPan+k,numPan+k) = -sum(L(kS(k),:) + L(kE(k),:)) + 2*pi;            % Vortex panel terms for each Kutta condition equation in A matrix
    end

    % Last element of b array (Kutta condition)
    for k = 1:1:numAF                                                           % Loop over all airfoils
        b(numPan+k) = -Vinf*2*pi*(sin(beta_Real(kS(k)))+sin(beta_Real(kE(k)))); % Set b array value from Kutta condition equation for each airfoil
    end

    % /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
    % ===== End/Enforce the Kutta Condition =====

    % Compute result array
    resArr = A\b;                                                               % Solution array including lambda and gamma

    % Separate lambda and gamma values from result array
    lambda = resArr(1:end-numAF);                                               % Source strengths for each panel
    gamma  = resArr(end-numAF+1:end);                                           % Vortex strength for each airfoil

    gammaArr = zeros(numPan,1);                                                 % Initialize the gamma array 
    for k = 1:1:numAF                                                           % Loop over every airfoil
        gammaArr(kS(k):kE(k)) = gamma(k);                                       % Set vortex strength for each panel on each airfoil
    end

    sumsumG = sum(gammaArr.*S_Real);                                            % Sum of all the vortex strengths (times panel lengths)

    % Chord length: Between minimum X coordinate and maximum X coordinate
    chord = max(XBT) - min(XBT);                                                % Difference between max X value and min X value
    Cl = 2*sumsumG/(chord*Vinf);                                                % Kutta-Joukowski lift Coefficient(chord length using chord)


    %% COMPUTE PANEL VELOCITIES AND PRESSURE COEFFICIENTS

    % Compute velocities on each panel
    Vt = zeros(numPan,1);                                                       % Initialize tangential velocity
    Cp = zeros(numPan,1);                                                       % Initialize pressure coefficient
    for i = 1:1:numPan
        term1 = Vinf*sin(beta_Real(i));                                         % Freestream term
        term2 = sum(lambda.*J(i,:)')/(2*pi);                                    % Source panel terms
        term3 = gammaArr(i)/2;                                                  % i == j vortex panel term
        term4 = -(gammaArr(i)/(2*pi)).*sum(L(i,:));                             % i ~= j vortex panel terms

        Vt(i) = term1 + term2 + term3 + term4;                                  % Compute total tangential velocity on panel i
        Cp(i) = 1-(Vt(i)/Vinf)^2;                                               % Compute pressure coefficient on panel i
    end

    
    %% Plot Data
    %% COMPUTE STREAMLINES - REFS [4,8]
        % Grid parameters
        xVals  = [min(XB)-0.5 max(XB)+0.5];                                     % X-grid extents [min, max]
        yVals  = [min(YB)-0.3 max(YB)+0.3];                                     % Y-grid extents [min, max]
        
        % Generate the grid points    
        Xgrid   = linspace(xVals(1),xVals(2),nGridX)';                          % X-values in evenly spaced grid        
        Ygrid   = linspace(yVals(1),yVals(2),nGridY)';                          % Y-values in evenly spaced grid 
        [XX,YY] = meshgrid(Xgrid,Ygrid);                                        % Create meshgrid from X and Y grid arrays

        % Initialize velocities
        Vx = zeros(nGridX,nGridY);                                              % Initialize X velocity matrix
        Vy = zeros(nGridX,nGridY);                                              % Initialize Y velocity matrix

        % Find which points are in/on the airfoil
        inon = zeros(nGridX,nGridY);                                            % Initialize in/on matrix with zeros
        for m = 1:1:nGridX                                                      % Loop over all X grid points
            for n = 1:1:nGridY                                                  % Loop over all Y grid points
                inMax = 0;                                                      % Set value to zero before checking all airfoils
                onMax = 0;                                                      % Set value to zero before checking all airfoils
                for k = 1:1:numAF                                               % Loop over all airfoils
                    [in,on] = inpolygon(XX(m,n),YY(m,n),XB_AF{k},YB_AF{k});     % Call function to determine if point is in/on the airfoil
                    inMax = max(inMax,in);                                      % Update value based on new information for current airfoil
                    onMax = max(onMax,on);                                      % Update value based on new information for current airfoil
                end
                if (inMax == 1 || onMax == 1)                                   % If the point is in/on any of the airfoils
                    inon(m,n) = 1;                                              % Set the value equal to one (otherwise it will stay zero)
                end
            end
        end

        % Solve for grid point X and Y velocities
        string = sprintf('%s\n\nComputing streamlines\n',string);               % Display status in Command Window
        set(handles.Dis1,'String',string);
        index = size(get(handles.Dis1,'string'), 1);                            % Get how many items are in the list box
        set(handles.Dis1,'Value',index);
        drawnow;
        
        for m = 1:1:nGridX                                                      % Loop over all X grid points
            index = size(get(handles.Dis1,'string'), 1);                        % Get how many items are in the list box
            set(handles.Dis1,'Value',index);
            string = sprintf('%s\n::Cal : X grid index: %i/%i\n',string,m,nGridX); 
            set(handles.Dis1,'String',string);                                  % Display current iteration in Command Window
            drawnow;
   
            for n = 1:1:nGridY                                                  % Loop over all Y grid points
                XP = XX(m,n);                                                   % Current iteration's X grid point
                YP = YY(m,n);                                                   % Current iteration's Y grid point
                [Mx,My] = STREAMLINE_SPM_N(XP,YP,XB,YB,phi,S,numPan,jInd);      % Compute Mx and My geometric integrals, REF [4]
                [Nx,Ny] = STREAMLINE_VPM_N(XP,YP,XB,YB,phi,S,numPan,jInd);      % Compute Nx and Ny geometric integrals, REF [8]

                if (inon(m,n) == 0)                                             % If the grid point is not in/on any airfoil
                    term1 = Vinf*cosd(AoA) + sum(lambda.*Mx/(2*pi));            % Uniform flow term and contribution from source panels
                    term2 = 0;                                                  % Reset value to zero for every grid point
                    for k = 1:1:numAF                                           % Loop over all airfoils
                        term2 = term2 + sum(-gamma(k).*Nx(kS(k):kE(k))/(2*pi)); % Contribution from vortex panels
                    end
                    Vx(m,n) = term1 + term2;                                    % Combine two terms to get X-direction velocity

                    term3 = Vinf*sind(AoA) + sum(lambda.*My/(2*pi));            % Uniform flow term and contribution from source panels
                    term4 = 0;                                                  % Reset value to zero for every grid point
                    for k = 1:1:numAF                                           % Loop over all airfoils
                        term4 = term4 + sum(-gamma(k).*Ny(kS(k):kE(k))/(2*pi)); % Contribution from source panels
                    end
                    Vy(m,n) = term3 + term4;                                    % Combine two terms to get Y-direction velocity
                end
            end
        end
    
        % Compute grid point velocity magnitude and pressure coefficient
        Vxy  = sqrt(Vx.^2 + Vy.^2);                                             % Compute magnitude of velocity vector []
        CpXY = 1-(Vxy./Vinf).^2;                                                % Pressure coefficient []
        
        string = sprintf('%s\n::Cal : Complete\nKutta-Joukowski lift Coefficient = %0.4f\n',...
                            string,Cl); 
        set(handles.Dis1,'String',string);                                      % Display current iteration in Command Window
        index = size(get(handles.Dis1,'string'), 1);                            % Get how many items are in the list box
        set(handles.Dis1,'Value',index);
%%
%Plot Data
%Define Plot type
Str = get(handles.popupmenu1 , 'String');
Val = get(handles.popupmenu1 , 'Value');
Plot_type = Str{Val};        

xVals  = [min(XB)-0.5 max(XB)+0.5];                                         % X-grid extents [min, max]
yVals  = [min(YB)-0.3 max(YB)+0.3];                                         % Y-grid extents [min, max]
    
% Streamline parameters       
stepsize = 0.01;                                                            % Step size for streamline propagation        
maxVert  = nGridX*nGridY*100;                                               % Maximum vertices        
slPct    = 50;                                                              % Percentage of streamlines of the grid        
Ysl      = linspace(yVals(1),yVals(2),floor((slPct/100)*nGridY))';          % Create array of Y streamline starting points
        

if strcmp(Plot_type,'Streamlines')
    cla(handles.Dis3,'reset'); hold(handles.Dis3,'on'); 
    set(handles.Dis3,'Color','White');                                                  % Set color to white
    set(handles.Dis3,'FontSize',12);                                                	% Set font size
    
    xlabel(handles.Dis3,'X Units');                                                     % Set X-label
    ylabel(handles.Dis3,'Y Units');                                                     % Set Y-label
    
    xlim(handles.Dis3,xVals);                                                           % Set axes equal,xVals);                                                           
    ylim(handles.Dis3,yVals);                                                           % Set axes equal,yVals);                                                            
    
    axis(handles.Dis3,'equal');                                                         % Set axes equal
    zoom(handles.Dis3,'reset');                                                         % Reset zoom 
    
    for i = 1:1:length(Ysl)                                                             % Loop over all Y streamline starting points
        streamline(handles.Dis3,XX,YY,Vx,Vy,xVals(1),Ysl(i),[stepsize,maxVert]);        % Plot the streamline
        set(handles.Dis3,'LineWidth',2);                                                % Set streamline line width
    end

    for i = 1:1:numAF                                                                   % Loop over all airfoils
        fill(handles.Dis3,XB_AF{i},YB_AF{i},'k');                                       % Plot airfoil bodies
    end  
    
elseif strcmp(Plot_type,'Pressure Contour')
    cla(handles.Dis3,'reset'); hold(handles.Dis3,'on'); 
    set(handles.Dis3,'Color','White');                                                  % Set color to white
    set(handles.Dis3,'FontSize',12);                                                	% Set font size
    
    xlabel(handles.Dis3,'X Units');                                                     % Set X-label
    ylabel(handles.Dis3,'Y Units');                                                     % Set Y-label
    
    xlim(handles.Dis3,xVals);                                                           % Set axes equal,xVals);                                                           
    ylim(handles.Dis3,yVals);                                                           % Set axes equal,yVals);                                                            
    
    axis(handles.Dis3,'equal');                                                         % Set axes equal
    zoom(handles.Dis3,'reset');                                                         % Reset zoom 
    
    contourf(handles.Dis3,XX,YY,CpXY,100,'EdgeColor','none');                           % Plot Cp contour
    for i = 1:1:numAF                                                                   % Loop over all airfoils
        fill(handles.Dis3,XB_AF{i},YB_AF{i},'k');                                       % Plot airfoil bodies
    end
    
elseif strcmp(Plot_type,'Pressure Coefficient Distribution')
    cla(handles.Dis3,'reset'); hold(handles.Dis3,'on'); 
    set(handles.Dis3,'Color','White');                                                  % Set color to white
    set(handles.Dis3,'FontSize',12);                                                	% Set font size
    
    set(handles.Dis3,'Ydir','reverse')                                                  % Reverse direction of Y-axis
    
    xlabel(handles.Dis3,'X/C');                                                         % Set X-label
    ylabel(handles.Dis3,'Cp');                                                          % Set Y-label                                                           
    
    zoom(handles.Dis3,'reset');                                                         % Reset zoom 
    
    XC = XC(panInd);
    plot(handles.Dis3,XC,Cp,'ks','MarkerFaceColor','k');                                % Plot Cp for all airfoils
    set(handles.Dis3,'UserData',get(handles.Dis3,'Position'))
end

%Show Selected Airfoil
XYB = zeros(max(numPtsAF),2,numAF);
    
for i=1:numAF        
    XYB(:,1,i) = cell2mat(XB_AF(i));        
    XYB(:,2,i) = cell2mat(YB_AF(i));        
    XYB(:,3,i) = 0;    
end
    
cla(handles.Dis2,'reset');                                                  % Get ready for plotting
set(handles.Dis2,'FontSize',12);                                            % Set font size
zoom(handles.Dis2,'reset');                                                 % Reset zoom

for i=1:numAF
    plot(handles.Dis2,XYB(:,1,i),XYB(:,2,i),'white');                       % Plot Airfoi
    hold(handles.Dis2,'on');
end 
set(handles.Dis2,'Color','black');                                          % Set color 

Information_string = sprintf('Airfoil : %s\n V_inf : %0.2f\n AOA : %0.2f\n nGrid : %d',char(AF_name),Vinf,AF_angle,nGridX);

set(handles.Information,'String',Information_string);
 

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Global Variable
global load_success cal_success
global numAF 
global xVals yVals stepsize maxVert Ysl
global XX YY XB_AF YB_AF XC
global Vx Vy Cp CpXY

if isempty(load_success)
    set(handles.Dis1,'Value',1);
    set(handles.Dis1,'String','::Error : Please "Run" First'); 
    return;
elseif load_success == 0   
    set(handles.Dis1,'Value',1);
    set(handles.Dis1,'String','::Error : Load Failure'); 
    return;    
elseif cal_success == 0     
    set(handles.Dis1,'Value',1);
    set(handles.Dis1,'String','::Error : Cal. Failure');
    return;    
end

%Update Plot type
Str = get(handles.popupmenu1 , 'String');
Val = get(handles.popupmenu1 , 'Value');
Plot_type = Str{Val};              

if strcmp(Plot_type,'Streamlines')
    cla(handles.Dis3,'reset'); hold(handles.Dis3,'on'); 
    set(handles.Dis3,'Color','White');                                                  % Set color to white
    set(handles.Dis3,'FontSize',12);                                                	% Set font size
    
    xlabel(handles.Dis3,'X Units');                                                     % Set X-label
    ylabel(handles.Dis3,'Y Units');                                                     % Set Y-label
    
    xlim(handles.Dis3,xVals);                                                           % Set axes equal,xVals);                                                           
    ylim(handles.Dis3,yVals);                                                           % Set axes equal,yVals);                                                            
    
    axis(handles.Dis3,'equal');                                                         % Set axes equal
    zoom(handles.Dis3,'reset');                                                         % Reset zoom 
    
    for i = 1:1:length(Ysl)                                                             % Loop over all Y streamline starting points
        streamline(handles.Dis3,XX,YY,Vx,Vy,xVals(1),Ysl(i),[stepsize,maxVert]);    % Plot the streamline
        set(handles.Dis3,'LineWidth',2);                                                % Set streamline line width
    end

    for i = 1:1:numAF                                                                   % Loop over all airfoils
        fill(handles.Dis3,XB_AF{i},YB_AF{i},'k');                                       % Plot airfoil bodies
    end  
    
elseif strcmp(Plot_type,'Pressure Contour')
    cla(handles.Dis3,'reset'); hold(handles.Dis3,'on'); 
    set(handles.Dis3,'Color','White');                                                  % Set color to white
    set(handles.Dis3,'FontSize',12);                                                	% Set font size
    
    xlabel(handles.Dis3,'X Units');                                                     % Set X-label
    ylabel(handles.Dis3,'Y Units');                                                     % Set Y-label
    
    xlim(handles.Dis3,xVals);                                                           % Set axes equal,xVals);                                                           
    ylim(handles.Dis3,yVals);                                                           % Set axes equal,yVals);                                                            
    
    axis(handles.Dis3,'equal');                                                         % Set axes equal
    zoom(handles.Dis3,'reset');                                                         % Reset zoom 
    
    contourf(handles.Dis3,XX,YY,CpXY,100,'EdgeColor','none');                           % Plot Cp contour
    for i = 1:1:numAF                                                                   % Loop over all airfoils
        fill(handles.Dis3,XB_AF{i},YB_AF{i},'k');                                       % Plot airfoil bodies
    end
    
elseif strcmp(Plot_type,'Pressure Coefficient Distribution')
    cla(handles.Dis3,'reset'); hold(handles.Dis3,'on'); 
    set(handles.Dis3,'Color','White');                                                  % Set color to white
    set(handles.Dis3,'FontSize',12);                                                	% Set font size
    
    set(handles.Dis3,'Ydir','reverse')                                                  % Reverse direction of Y-axis
    
    xlabel(handles.Dis3,'X/C');                                                         % Set X-label
    ylabel(handles.Dis3,'Cp');                                                          % Set Y-label                                                           
    
    zoom(handles.Dis3,'reset');                                                         % Reset zoom 
    
    plot(handles.Dis3,XC,Cp,'ks','MarkerFaceColor','k');                                % Plot Cp for all airfoils
    set(handles.Dis3,'UserData',get(handles.Dis3,'Position'))
end



function Information_Callback(hObject, eventdata, handles)
% hObject    handle to Information (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Information as text
%        str2double(get(hObject,'String')) returns contents of Information as a double


% --- Executes during object creation, after setting all properties.
function Information_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Information (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
