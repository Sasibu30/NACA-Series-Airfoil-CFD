function [XB,YB,numPtsAF,numPanAF,title,string,success] = check_input_gui(AF_name,AF_load,AF_flip,AF_scale,AF_angle,AF_offset,Vinf)
   
    string = sprintf('============\n%s\n============\n','Start Program');   
    % Number of airfoils    
    numAF = length(AF_name);                                                    % Number of airfoils
    
    XB       = cell(numAF,1);                                                   % Initialize X boundary points
    YB       = cell(numAF,1);                                                   % Initialize Y boundary points
    numPtsAF = zeros(numAF,1);                                                  % Number of points
    numPanAF = zeros(numAF,1);                                                  % Number of panels
    
    string = sprintf('%s\n%s\n',string,'Validation Start');
    
    % Check number of airfoils
    if numAF == 1
        title = char('single_');
    elseif numAF == 2
        title = char('double_');
    elseif numAF == 3
        title = char('triple_');    
    else
        string = sprintf('%s\n::Error : We can not get more than 3 Airfoils Dat.\n',string);
        success = 0;
        return;
    end
    string = sprintf('%s\n::Msg : Number of Airfoils is [%d]',string,numAF);
    
    for i=1:numAF       
        title = char(strcat(title,AF_name(i),'_'));    
    end
    string = sprintf('%s\n::Msg : Title is "%s"',string,title);
    
    
    % Check Airfoil Get. Dat
    if ~exist('./Airfoil_DAT','dir')
        mkdir('./Airfoil_DAT');
    end
    for i = 1:1:numAF
        %Call Airfoil Dat.( The user should load a Airfoil DAT File )
        if (AF_load(i) == 0)                                                    % If do not use the NACA Series airfoill     
            airfoilName = char(strcat('./Airfoil_DAT/', AF_name(i), '.dat'));                                       
        elseif (AF_load(i) == 1)                                                % Set the NACA Series airfoill
            airfoilName = char(strcat('./Airfoil_DAT/naca', AF_name(i), '.dat'));
        end       
        fidAirfoil = fopen(airfoilName,'r');
        
        %Check Success of Dat. File Loading 
        if fidAirfoil > 0
            success = 1;
        else
            string = sprintf('%s\n::Error : Dat. Load Failure',string);
            string = sprintf('%s\n::Error : Check Name of Airfoil',string);
            string = sprintf('%s\n::Error : or',string);
            string = sprintf('%s\n::Error : Check "%s"',string,airfoilName);
            string = sprintf('%s\n::Error : Critical Error about Load Dat',string);
            success = 0;
            return;
        end

        %Load Data
        dataBuffer = textscan(fidAirfoil,'%f %f',...
                'CollectOutput',1,'Delimiter','','HeaderLines',1);              % Load selected airfoil

        xFoilResults.XB = dataBuffer{1}(:,1);                                   % Airfoil boundary X-points
        xFoilResults.YB = dataBuffer{1}(:,2);                                   % Airfoil boundary Y-points

        XB{i,1} = AF_scale(i)*xFoilResults.XB;                                  % Boundary point X-coordinate
        YB{i,1} = AF_flip(i)*AF_scale(i)*xFoilResults.YB;                       % Boundary point Y-coordinate

        fclose(fidAirfoil);                                                     % Close the airfoil file
        
        % Rotate the airfoil
        v      = [XB{i} YB{i}];                                                 % Concatenate for rotation
        R      = [cosd(-AF_angle(i)) -sind(-AF_angle(i));                       % Rotation matrix based on each airfoil's rotation angle (AF_angle)
                  sind(-AF_angle(i)) cosd(-AF_angle(i))];
        vo     = (R*v')';                                                       % Rotate the airfoil data
        x_rot  = vo(:,1);                                                       % Extract rotated X boundary points
        y_rot  = vo(:,2);                                                       % Extract rotated Y boundary points

        XB{i,1} = x_rot + AF_offset(i,1);                                       % Offset the X boundary points by user-defined placement (AF_offset)
        YB{i,1} = y_rot + AF_offset(i,2);                                       % Offset the Y boundary points by user-defined placement (AF_offset)

        numPtsAF(i,1) = length(XB{i});                                          % Number of points for each airfoil
        numPanAF(i,1) = length(XB{i})-1;                                        % Number of panels for each airfoil
    end   
    
    % Check Airfoil AOA
    if AF_angle > 10 || AF_angle < -10
        answer = questdlg('[ AOA > 10 or AOA < -10 ] can cause too much difference between actual and calculated data, Do you want to continue?',...
            'Confirm',...
            'Yes','No','No');

        % Handle response
        switch answer
            case 'Yes'
                success = 1;
            case 'No'
                success = -1;
                string = sprintf('::Msg : Program "Run" Cancel\n');
                return;
        end
    end
end