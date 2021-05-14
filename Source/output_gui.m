function[success] = output_gui(Plot_type,XB_AF,YB_AF,XC,Cp,XX,YY,Vx,Vy,CpXY,AF_name,numPtsAF,numAF,title)   
    if ~exist('./Output','dir')
        mkdir('./Output');
    end
    
    if strcmp(Plot_type,'Airfoil_Dat')
        saveFlnm    = ['./Output/Save_', title, 'Airfoil.dat'];                     % Airfoil coordinates save-to file
        XYB = zeros(max(numPtsAF),2,numAF);
        
        % Delete files if they exist
        if (exist(saveFlnm,'file'))                                                 % If airfoil coordinate file exists
            delete(saveFlnm);                                                       % Delete it
        end

        fidAirfoil1 = fopen(saveFlnm,'w');                                          % Open the airfoil file
        
        %check
        [a,b] = size(XYB);
        [A,B] = size(cell2mat(XB_AF));
        
        if fidAirfoil1 < 0 
            success = 0;
            fclose(fidAirfoil1);
            return     
        elseif a ~= A
            success = 0;
            fclose(fidAirfoil1);
            return     
        else
            success = 1;
        end
        
        for i=1:numAF
            XYB(:,1,i) = cell2mat(XB_AF(i));
            XYB(:,2,i) = cell2mat(YB_AF(i));
            XYB(:,3,i) = 0;
        end
       
         % Save the airfoil coordinate file
        for i=1:numAF
            fprintf(fidAirfoil1,'Zone T = "%s Airfoil Coordinate" \n',cell2mat(AF_name(i)));
            fprintf(fidAirfoil1,'%f, %f \n',XYB(:,1:2,i)');  
        end
        fclose(fidAirfoil1);
    elseif strcmp(Plot_type,'Pressure Contour')
        saveFlnmCp  = ['./Output/Save_', title, 'Pressure Contour.dat'];            % Pressure Contour save-to file

        % Delete files if they exist
        if (exist(saveFlnmCp,'file'))                                               % If airfoil Cp file exists
            delete(saveFlnmCp);                                                     % Delete it
        end

        fidAirfoil2 = fopen(saveFlnmCp,'w');                                        % Open the Cp file
        
        if fidAirfoil2 > 0
            success = 1;
        else
            fclose(fidAirfoil2);
            success = 0;
            return
        end
        
        [I,J] = size(XX); K=1; CpXXYY=zeros(I*J,3); inde=1;
        for i=1:I
            for j=1:J
                CpXXYY(inde,1) = XX(1,j)';
                CpXXYY(inde,2) = YY(i,1)';
                CpXXYY(inde,3) = CpXY(i,j);
                inde = inde+1;
            end
        end

        % Save the airfoil Cp file
        fprintf(fidAirfoil2,'VARIABLES="X" "Y" "Cp" \n');
        fprintf(fidAirfoil2,'Zone T = "%s Calculation Area Cp" \n',title);
        fprintf(fidAirfoil2,'I=%d, J=%d, K=%d F=point\n',I,J,K);
        fprintf(fidAirfoil2,'%f, %f, %f \r\n',CpXXYY');                             

        fclose(fidAirfoil2);
    elseif strcmp(Plot_type,'Streamlines')
        saveFlnmVec  = ['./Output/Save_', title, '_Streamlines.dat'];                       % Airfoil Vec save-to file

        % Delete files if they exist
        if (exist(saveFlnmVec,'file'))                                              % If airfoil Vec file exists
            delete(saveFlnmVec);                                                    % Delete it
        end

        fidAirfoil3 = fopen(saveFlnmVec,'w');                                       % Open the Vec file
        
        if fidAirfoil3 > 0
            success = 1;
        else
            success = 0;
            fclose(fidAirfoil3);
            return
        end
        
        [I,J] = size(CpXY);
        VecXXYY=zeros(I*J,4); inde=1;
        for i=1:I
            for j=1:J
                VecXXYY(inde,1) = XX(1,j)';
                VecXXYY(inde,2) = YY(i,1)';
                VecXXYY(inde,3) = Vx(i,j);
                VecXXYY(inde,4) = Vy(i,j);
                inde = inde+1;
            end
        end

        %save the airfoil Vec file
        fprintf(fidAirfoil3,'VARIABLES="X" "Y" "Vx" "Vy" \n');
        fprintf(fidAirfoil3,'Zone T = "%s Streamlines" \n',title);
        fprintf(fidAirfoil3,'I=%d, J=%d, K=%d \n',max(numPtsAF),1,1);
        fprintf(fidAirfoil3,'%f, %f, %f, %f \r\n',VecXXYY');

        fclose(fidAirfoil3);
        
    elseif strcmp(Plot_type,'Pressure Coefficient Distribution')
        saveFlnmCpDis  = ['./Output/Save_', title, 'Pressure Coefficient Distribution.dat'];           % Airfoil Cp Distribution save-to file
        
        % Delete files if they exist
        if (exist(saveFlnmCpDis,'file'))                                              % If Cp Distribution file exists
            delete(saveFlnmCpDis);                                                   % Delete it
        end

        fidAirfoil4 = fopen(saveFlnmCpDis,'w');                                       % Open Cp Distribution file
        
        if fidAirfoil4 > 0
            success = 1;
        else
            fclose(fidAirfoil4);
            success = 0;
            return
        end
        
        XC_CP(:,1) = XC;
        XC_CP(:,2) = Cp;
        fprintf(fidAirfoil4,'Zone T = "%s Cp Distribution at Airfoil" \n',title);
        fprintf(fidAirfoil4,'%f, %f \r\n',XC_CP');
        
        fclose(fidAirfoil4);
    end
 
end