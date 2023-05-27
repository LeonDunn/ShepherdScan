function [LogDate, ID, Plan, Beam, VarianSerial, pp_init, Result, RMS_MLC, MAX_X1, MAX_X2, MAX_Y1, MAX_Y2, MAX_LAT, MAX_LONG, MAX_VERT, MAX_COLL, MAX_GANTRY, fname_A, fname_E, fname_G, fname_out, BEAM_ON] = eFluenceQA_Tlog(app, filename, OutFolder, GammaDD, GammaDTA, GammaThresh, PassingCriteria)
%     close all
%     clear all
%     OutFolder = 'C:\Temp\Output\'
%     filename = '129934_N19a R Br_LN_1_1.5_20230104173137.bin'
%     GammaDD = 0.03
%     GammaDTA = 3
%     GammaThresh = 0.1
%     PassingCriteria = 0.95


    [folder, fname] = fileparts(filename);
    fname_out = fname;

    
    fid=fopen(filename);
    %
    int_skip = 20;
    % File specification from Varian - header is fixed 1024 bytes long
    % Groups\Gainey\MATLAB\TrueBeam 1.5 Trajectory Log File Specification.pdf
    
    % see tables on pages 7-13 for more details
       
    signature=char(fread(fid,16));
    
    %fprintf('\nEnd of signature found at position %d\n',ftell(fid) )
    
    version=char(fread(fid,16));
    
    %fprintf('\nEnd of version found at position %d\n',ftell(fid) )
    %headerSize=uint32(fread(fid,4,'integer*4')); % integer fixed 1024 bytes
    
     headerSize=uint32(fread(fid,1,'int32')); % integer fixed 1024 bytes
     samplingInterval=uint32(fread(fid,1,'int32')); % integer milliseconds
     numberOfAxesSampled=uint32(fread(fid,1,'int32')); % integer
     
    
    % Axis enumeration - table page 8
    
    axisEnumeration=uint32(fread(fid,numberOfAxesSampled,'int32'));
    
    %samplesPerAxis
    
    samplesPerAxis=uint32(fread(fid,numberOfAxesSampled,'int32'));
    axisScale=uint32(fread(fid,1,'int32'));
    
    %numberOfSubbeams
    numberOfSubBeams = (fread(fid,1,'int32'));
    
    % isTruncated
    isTruncated=uint32(fread(fid,1,'int32'));
    
    %numberOfSnapshots
    numberOfSnapShots=uint32(fread(fid,1,'int32'));
    
    %MLC model
    mlcModel=uint32(fread(fid,1,'int32'));
    if (mlcModel==2) 
        MLC='NDS 120';
    else 
        if (mlcModel==3) 
            MLC='NDS 120 HD';
        end
    end
    
    
    %reserved - Varian specific data stored here, lets skip over it to start of
    %subbeam data (position 1024)
    
    
    % now skip over first 1024 bytes of header: second argument is offset value
        
    fseek(fid,headerSize-1,'bof'); % move pointer 1024 from beginning of file 'bof' http://www.mathworks.de/de/help/matlab/ref/fseek.html
    
    subBeam(numberOfSubBeams)=struct('cp',[],'mu',[],'radTime',[],'seq',[],'name',[]); 
    % now read in all subBeams- each subBeam is currently fixed to 560 bytes
    % (Varian specification)
    
    for i=1:numberOfSubBeams
       %fprintf('\nStart of subBeam %d, position= %d\n',i,ftell(fid))
       
       subBeam(i).cp=(fread(fid,1,'*int32'));
       
       %fprintf('subBeam %d cp=%d\n',i,subBeam(i).cp)
    
       subBeam(i).mu=(fread(fid,1,'*single'));
       
       subBeam(i).radTime =fread(fid,1,'*single');
        
       %temp3=fread(fid,4,'uint8');
       %subBeam(i).seq=typecast(uint8(temp3),'int32');
       
       subBeam(i).seq=(fread(fid,1,'*int32'));
       
       subBeam(i).name = fread(fid,32,"char");
       %fprintf('\nEnd of .name position in subBeam %d, = %d\n',i,ftell(fid));
       
       %skip over Varian reserved data (32 bytes) at end of subBeam structure
       
       %varianReserve=char(fread(fid,32));
       % this line needs to edited for version 2.0 (bytes shift = 560bytes) not
       % 80 bytes
       %fseek(fid,1024+(i*80),'bof'); % skip over reserved 32 bytes i.e. count 80 bytes from end of header (1024) of current subBeam 1024+i*80
       if str2double(version')>=2.0
            fseek(fid,1024+(i*560),'bof'); %uncommment this line for version 2.0
       else
            fseek(fid,1024+(i*80),'bof');
       end
       %i.e. skip 560 bytes per subBeam from end of header (1024)
    end
    
    % snapShot structure consists of paired values: expected (E) and actual (A) as
    % described in Varian specification
    % axisEnumeration yields the following 
    % [0;1;2;3;4;5;6;7;8;9;40;41;42;50] corresponding to
    
    % 0 - Coll Rtn
    % 1 - Gantry Rtn
    % 2 - Y1
    % 3 - Y2
    % 4 - X1
    % 5 - X2
    % 6 - Couch Vrt
    % 7 - Couch Lng
    % 8 - Couch Lat
    % 9 - Couch Rtn
    % 40 - MU
    % 41 - Beam Hold
    % 42 - Control Point
    % 50 - MLC
    
    Leaves=samplesPerAxis(numberOfAxesSampled)-2;
    
    snapShot(numberOfSnapShots)=struct('colRotationE',[],'colRotationA',[], ...
                   'gantryRotationE',[],'gantryRotationA',[],...
                    'Y1_E',[],'Y1_A',[],'Y2_E',[],'Y2_A',[],...
                    'X1_E',[],'X1_A',[],'X2_E',[],'X2_A',[],...
                    'couchVrtE', [], 'couchVrtA', [],...
                    'couchLngE', [], 'couchLngA', [],...
                    'couchLatE',[], 'couchLatA',[],...
                    'couchRotationE', [], 'couchRotationA', [],...
                    'MU_E', [],'MU_A', [], ...
                    'beamHoldE', [], 'beamHoldA', [], ...
                    'controlPointE', [], 'controlPointA', [], ...
                    'CarrA_E', [], 'CarrA_A', [], ...
                    'CarrB_E', [], 'CarrB_A', [], ...
                    'MLC_E', {zeros(1,Leaves)}, 'MLC_A', {zeros(1,Leaves)});
    
    %Now comes the hard bit, reading in all the snapshots for each axis
    %(samplesPerAxis), Actual and Expected values
    
    totalNumberOfSamples=sum(samplesPerAxis); % summation over all axes
    
    % first read in all float arrays, and then split into respective snapShots
    % initialise temp unstructured array to store data
    
    % fieldNames={'colRotationE','colRotationA', ...
    %                'gantryRotationE','gantryRotationA',...
    %                 'Y1_E','Y1_A',...
    %                 'Y2_E','Y2_A',...
    %                 'X1_E','X1_A',...
    %                 'X2_E','X2_A',...
    %                 'couchVrtE', 'couchVrtA',...
    %                 'couchLngE', 'couchLngA',...
    %                 'couchLatE', 'couchLatA',...
    %                 'couchRotationE', 'couchRotationA',...
    %                 'MU_E', 'MU_A',...
    %                 'beamHoldE', 'beamHoldA',...
    %                 'controlPointE','controlPointA', ...
    %                 'MLC_E', 'MLC_A' };
                
    % here we just grab all axis data for current (i) snapshot into temporary
    % cellNum array temp5. After we have all information we can assign data to a
    % structured array (snapShot)
    
    %temp6=zeros(1,2*totalNumberOfSamples);
    
    % temp5={zeros(1,2*totalNumberOfSamples)};
     
    beamOn=zeros(numberOfSnapShots,1);
    
    
    %lastCPArc(numberOfSubBeams)=subBeam(numberOfSubBeams-1).cp+subBeam(numberOfSubBeams).cp;
    
    
    % b(numberOfSubBeams)=struct('On',[],'Off',[]);
    beamCP(numberOfSubBeams)=struct('On',[],'Off',[]);
    
    j=1;
    
    for i=1:numberOfSnapShots
        % actual_(1,j,a1),expected_(1,j,a1),actual_(1,j,a2),expected_(1,j,a2),....
        % see figure on page 11 of Specification
    
          %temp5{i}= single(fread(fid,2*totalNumberOfSamples,'single'));
    
          temp5{i} = single(fread(fid,2*totalNumberOfSamples,'*single'));
          %temp5{i} = double(temp5{i});
          BH = int32(temp5{i}(28));
          
          %if the beam is held, we just want to skip recording the data and
          %move on
          if (BH >= 1) 
              %disp('beam holding')
              temp5{i} = [];
              continue
          end
          
         
    % need to convert angles (colli, gantry, couch) from intrinsic Varian scale (+1°) to IEC 61217 scale (+179°)
     
    % snapShot(i).colRotationE=single(temp5{i}(1));
          snapShot(j).colRotationE = Varian2IEC((temp5{i}(1)));
          
          snapShot(j).colRotationA = Varian2IEC((temp5{i}(2)));
        
          snapShot(j).gantryRotationE = Varian2IEC((temp5{i}(3)));
          
          snapShot(j).gantryRotationA = Varian2IEC((temp5{i}(4)));
          
          snapShot(j).Y1_E = (temp5{i}(5));
          
          snapShot(j).Y1_A = (temp5{i}(6));
          
          snapShot(j).Y2_E = (temp5{i}(7));
          
          snapShot(j).Y2_A = (temp5{i}(8));
          
          snapShot(j).X1_E = (temp5{i}(9));
          
          snapShot(j).X1_A = (temp5{i}(10));
          
          snapShot(j).X2_E = (temp5{i}(11));
          
          snapShot(j).X2_A = (temp5{i}(12));
          
          snapShot(j).couchVrtE = (temp5{i}(13));
          
          snapShot(j).couchVrtA = (temp5{i}(14));
          
          snapShot(j).couchLngE = (temp5{i}(15));
          
          snapShot(j).couchLngA = (temp5{i}(16));
          
          snapShot(j).couchLatE = (temp5{i}(17));
          
          snapShot(j).couchLatA = (temp5{i}(18));
          
          snapShot(j).couchRotationE = Varian2IEC((temp5{i}(19)));
          
          snapShot(j).couchRotationA = Varian2IEC((temp5{i}(20)));
    
          snapShot(j).MU_E = (temp5{i}(25));
          
          snapShot(j).MU_A = (temp5{i}(26));
          
          snapShot(j).beamHoldE = (temp5{i}(27));
          
          snapShot(j).beamHoldA = (temp5{i}(28));
          
          snapShot(j).controlPointE = (temp5{i}(29));
          
          snapShot(j).controlPointA = (temp5{i}(30));
          
          snapShot(j).CarrA_E = (temp5{i}(31));
          
          snapShot(j).CarrA_A = (temp5{i}(32));
          
          snapShot(j).CarrB_E = (temp5{i}(33));
          
          snapShot(j).CarrB_A = (temp5{i}(34));
          
          %last butnot least comes the MLC
          
          snapShot(j).MLC_E = (temp5{i}(35:2:end-1));
          snapShot(j).MLC_A = (temp5{i}(36:2:end));
    
          j = j+1;
    end
    
    %remove all the beams with zero MU
    
   snapShotFinal = snapShot(all(~cellfun(@isempty,struct2cell(snapShot))));
   
    



    headerInfo = struct('signature', signature, ...
                'version', version, ...
                'headerSize', headerSize, ...
                'samplingInterval', samplingInterval, ...
                'numberOfAxesSampled',numberOfAxesSampled,...
                'axisEnumeration',axisEnumeration, ...
                'samplesPerAxis',samplesPerAxis,...
                'axisScale', axisScale, ...
                'numberOfSubBeams', numberOfSubBeams, ...
                'isTruncated',isTruncated, ...
                'numberOfSnapShots', numberOfSnapShots, ...
                'MLC', MLC,'totalSamples', totalNumberOfSamples );
      
    
    fclose(fid);
    
    %1. Check for subbeams
    size_subBeam = size(subBeam,2);
    
    if size_subBeam > 1 %sub beams exist
    i = 0;
      
        %get the value for the cp which defines the starting poing of the beam
        CP_ind = [subBeam.cp];
        CP_ind = single(CP_ind);
        
    CP_list = unique([snapShotFinal.controlPointA]);
        
    % Give "true" if the element in "a" is a member of "b".
    [~, c] = ismember(CP_list, CP_ind);
    % Extract the elements of a at those indexes.
    indexes = find(c,size_subBeam+1,'last');
    size_indexes = size(indexes,2);


    if size_indexes > size_subBeam
        diff_size = size_indexes - size_subBeam;
        indexes(end-diff_size+1:end) = [];
    end
    
        for j = 1:size(snapShotFinal,2) %sample
                    
                    ControlPoint(j,1) = snapShotFinal(j).controlPointA;
                    
                    if  ismember(j,indexes)
                        
                            i = i+1;
                    end
                    
                    %actual MLC               
                    CollAngle(j,i) = snapShotFinal(j).colRotationA;
                    
                    X1_JAW{i}{j} = snapShotFinal(j).X1_A*10;
                    X2_JAW{i}{j} = snapShotFinal(j).X2_A*10;
                    Y1_JAW{i}{j} = snapShotFinal(j).Y1_A*10;
                    Y2_JAW{i}{j} = snapShotFinal(j).Y2_A*10;
                   
                    MLC_A{i}{j} = flipud(snapShotFinal(j).MLC_A(1:60,1)*10);
                    MLC_B{i}{j} = flipud(snapShotFinal(j).MLC_A(61:end,1)*10);
                    MU_seg(j,i) = snapShotFinal(j).MU_A;

                    CollAnglee(j,i) = snapShotFinal(j).colRotationE;
                    
                    X1_JAWe{i}{j} = snapShotFinal(j).X1_E*10;
                    X2_JAWe{i}{j} = snapShotFinal(j).X2_E*10;
                    Y1_JAWe{i}{j} = snapShotFinal(j).Y1_E*10;
                    Y2_JAWe{i}{j} = snapShotFinal(j).Y2_E*10;
                   
                    MLC_Ae{i}{j} = flipud(snapShotFinal(j).MLC_E(1:60,1)*10);
                    MLC_Be{i}{j} = flipud(snapShotFinal(j).MLC_E(61:end,1)*10);
                    MU_sege(j,i) = snapShotFinal(j).MU_E;


                



        end
        
    else
    i = 1;

        for j = 1:size(snapShotFinal,2) %sample
                       
                    ControlPoint(j,i) = snapShotFinal(j).controlPointA;
                    CollAngle(j,i) = snapShotFinal(j).colRotationA;
                 
                    X1_JAW{i}{j} = snapShotFinal(j).X1_A*10;
                    X2_JAW{i}{j} = snapShotFinal(j).X2_A*10;
                    Y1_JAW{i}{j} = snapShotFinal(j).Y1_A*10;
                    Y2_JAW{i}{j} = snapShotFinal(j).Y2_A*10;
                   
                    MLC_A{i}{j} = flipud(snapShotFinal(j).MLC_A(1:60,1)*10);
                    MLC_B{i}{j} = flipud(snapShotFinal(j).MLC_A(61:end,1)*10);
                    MU_seg(j,i) = snapShotFinal(j).MU_A;
                    
                    %expected
                    CollAnglee(j,i) = snapShotFinal(j).colRotationE;
                    
                    X1_JAWe{i}{j} = snapShotFinal(j).X1_E*10;
                    X2_JAWe{i}{j} = snapShotFinal(j).X2_E*10;
                    Y1_JAWe{i}{j} = snapShotFinal(j).Y1_E*10;
                    Y2_JAWe{i}{j} = snapShotFinal(j).Y2_E*10;
                   
                    MLC_Ae{i}{j} = flipud(snapShotFinal(j).MLC_E(1:60,1)*10);
                    MLC_Be{i}{j} = flipud(snapShotFinal(j).MLC_E(61:end,1)*10);
                    MU_sege(j,i) = snapShotFinal(j).MU_E;
                    
        end
    
    end
    
    

    
    %minutes
    %decimate by 5's
    [rows cols] = size(MLC_A)
    for i = 1:cols
        BeamOnTime(i,1) = duration(seconds((length(MLC_A{i})*20)/1000))
    end

    MU_seg = MU_seg(1:int_skip:end,:);
    MU_sege = MU_sege(1:int_skip:end,:);
    CollAnglee = CollAnglee(1:int_skip:end,:);
    CollAngle = CollAngle(1:int_skip:end,:);
        
    if size_subBeam == 1
        MLC_A{:} = MLC_A{:}(1:int_skip:end);
        MLC_Ae{:} = MLC_Ae{:}(1:int_skip:end);
        MLC_B{:} = MLC_B{:}(1:int_skip:end);
        MLC_Be{:} = MLC_Be{:}(1:int_skip:end);
        X1_JAW{:} = X1_JAW{:}(1:int_skip:end);
        X2_JAW{:} = X2_JAW{:}(1:int_skip:end);
        Y1_JAW{:} = Y1_JAW{:}(1:int_skip:end);
        Y2_JAW{:} = Y2_JAW{:}(1:int_skip:end);
        X1_JAWe{:} = X1_JAWe{:}(1:int_skip:end);
        X2_JAWe{:} = X2_JAWe{:}(1:int_skip:end);
        Y1_JAWe{:} = Y1_JAWe{:}(1:int_skip:end);
        Y2_JAWe{:} = Y2_JAWe{:}(1:int_skip:end);
    end
    
    


    if size_subBeam > 1
        for i = 1:size_subBeam
            
           

            if isempty(MLC_A{i})
                continue
            end

            MLC_A{i} = MLC_A{i}(1:int_skip:end);
            MLC_Ae{i} = MLC_Ae{i}(1:int_skip:end);
            MLC_B{i} = MLC_B{i}(1:int_skip:end);
            MLC_Be{i} = MLC_Be{i}(1:int_skip:end);
            X1_JAW{i} = X1_JAW{i}(1:int_skip:end);
            X2_JAW{i} = X2_JAW{i}(1:int_skip:end);
            Y1_JAW{i} = Y1_JAW{i}(1:int_skip:end);
            Y2_JAW{i} = Y2_JAW{i}(1:int_skip:end);
            X1_JAWe{i} = X1_JAWe{i}(1:int_skip:end);
            X2_JAWe{i} = X2_JAWe{i}(1:int_skip:end);
            Y1_JAWe{i} = Y1_JAWe{i}(1:int_skip:end);
            Y2_JAWe{i} = Y2_JAWe{i}(1:int_skip:end);
        end
    end

 
    
    %correct the coll angles that are counter clockwise
    MU_seg = [repmat(0,[1,size(MU_seg,2)]); diff(MU_seg)];
    MU_seg(MU_seg<0) = 0;
    MU_seg(MU_seg>100) = 0;

    MU_sege = [repmat(0,[1,size(MU_sege,2)]); diff(MU_sege)];
    MU_sege(MU_sege<0) = 0;
    MU_sege(MU_sege>100) = 0;

    fluence_grid = ones(400,400);
    
        
    if mlcModel == 2 %MMLC
            n_large_leaves = 10;
            n_small_leaves = 40;
            leaf_width_large = 10; %mm
            leaf_width_small = 5; %mm
            large_leaf_offset = leaf_width_large;
            leaf_length = 200;
            num_leaves = 60;
            
               %no graphics
                j = 1;
                %create the rectangle MLC leaves
                for i = 0:10:90
                    A(j) = images.roi.Rectangle('Position',[0,i,200,leaf_width_large],'Color','r','InteractionsAllowed','none');
                    A(j).FaceAlpha = 0.9;
                    j=j+1;
                end
                %create the rectangle MLC leaves
                for i = 100:5:295
                    A(j) = images.roi.Rectangle('Position',[0,i,200,leaf_width_small],'Color','r','InteractionsAllowed','none');
                    j=j+1;
                end
                %create the rectangle MLC leaves
                for i = 300:10:390
                    A(j) = images.roi.Rectangle('Position',[0,i,200,leaf_width_large],'Color','r','InteractionsAllowed','none');
                    A(j).FaceAlpha = 0.9;
                    j=j+1;
                end
                    
                j=1;
                %create the rectangle MLC leaves
                for i = 0:10:90
                    B(j) = images.roi.Rectangle('Position',[200,i,200,leaf_width_large],'Color','b','InteractionsAllowed','none');
                    B(j).FaceAlpha = 0.9;
                    j=j+1;
                end
                %create the rectangle MLC leaves
                for i = 100:5:295
                    B(j) = images.roi.Rectangle('Position',[200,i,200,leaf_width_small],'Color','b','InteractionsAllowed','none');
                    B(j).FaceAlpha = 0.9;
                    j=j+1;
                end
                %create the rectangle MLC leaves
                for i = 300:10:390
                    B(j) = images.roi.Rectangle('Position',[200,i,200,leaf_width_large],'Color','b','InteractionsAllowed','none');
                    B(j).FaceAlpha = 0.9;
                    j=j+1;
                end
    end
    if mlcModel == 3 %HDMLC

            n_large_leaves = 14;
            n_small_leaves = 60;
            leaf_width_large = 5; %mm
            leaf_width_small = 2.5; %mm
            large_leaf_offset = leaf_width_large;
            leaf_length = 200;
            num_leaves = 60;
           
            centre_offset = 82.75; %this is the shift to make the MLC centred on the plane
             %no graphics
            j = 1;
            %create the rectangle MLC leaves
            for i = 0+centre_offset:5:70+centre_offset
                A(j) = images.roi.Rectangle('Position',[0,i,200,leaf_width_large],'Color','r','InteractionsAllowed','none');
                A(j).FaceAlpha = 0.1;
                j=j+1;
            end
            %create the rectangle MLC leaves
            for i = 75+centre_offset:2.5:155+centre_offset
                A(j) = images.roi.Rectangle('Position',[0,i,200,leaf_width_small],'Color','r','InteractionsAllowed','none');
                A(j).FaceAlpha = 0.1;
                j=j+1;
                
            end
            %create the rectangle MLC leaves
            for i = 157.5+centre_offset:5:227.5+centre_offset
        
                A(j) = images.roi.Rectangle('Position',[0,i,200,leaf_width_large],'Color','r','InteractionsAllowed','none');
                A(j).FaceAlpha = 0.1;
                j=j+1;
            end
                
            j=1;
            %create the rectangle MLC leaves
            for i = 0+centre_offset:5:70+centre_offset
                B(j) = images.roi.Rectangle('Position',[200,i,200,leaf_width_large],'Color','b','InteractionsAllowed','none');
                B(j).FaceAlpha = 0.1;
                j=j+1;
            end
            %create the rectangle MLC leaves
            for i = 75+centre_offset:2.5:155+centre_offset
                B(j) = images.roi.Rectangle('Position',[200,i,200,leaf_width_small],'Color','b','InteractionsAllowed','none');
                B(j).FaceAlpha = 0.1;
                j=j+1;
            end
            %create the rectangle MLC leaves
            for i =  157.5+centre_offset:5:227.5+centre_offset
                B(j) = images.roi.Rectangle('Position',[200,i,200,leaf_width_large],'Color','b','InteractionsAllowed','none');
                B(j).FaceAlpha = 0.1;
                j=j+1;
            end
    end
        
        %create the JAW ROI
        X1_JAW_ROI(1) = images.roi.Rectangle('Position',[0,0,200,400],'Color','c','InteractionsAllowed','none');
        X2_JAW_ROI(1) = images.roi.Rectangle('Position',[200,0,200,400],'Color','m','InteractionsAllowed','none');
        Y1_JAW_ROI(1) = images.roi.Rectangle('Position',[0,200,400,200],'Color','y','InteractionsAllowed','none');
        Y2_JAW_ROI(1) = images.roi.Rectangle('Position',[0,0,400,200],'Color','k','InteractionsAllowed','none');
       


        %set the original positions to reset to
        for i = 1:length(A)
            MLC_A_shift(i) = A(i).Position(1);
            MLC_B_shift(i) = B(i).Position(1);
        end
        
        %append the original position to a new var
        JAW_X1_shift(1) = X1_JAW_ROI(1).Position(1);
        JAW_X2_shift(1) = X2_JAW_ROI(1).Position(1);
        JAW_Y1_shift(1) = Y1_JAW_ROI(1).Position(2);
        JAW_Y2_shift(1) = Y2_JAW_ROI(1).Position(2);

        
        fluence_grid_total = zeros(400,400);
        MLC_mask = ones(400,400);
        final_fluence_plan_A = zeros(400,400);
    
        for k = 1:length(MLC_A)
        
            
                VarianSerial{k,1} = cell2mat(regexp(folder,'H\d*','match'));
            if isempty(VarianSerial{k})
                VarianSerial{k,1} = 'H000000'; %unknown machine
            end
 
            if isempty(MLC_A{k})
                continue
            end

         
            for l = 1:length(MLC_A{k})
                if isempty(MLC_A{k}{l})
                    continue
                end
                  
                  opening_pos{k}(:,1) = MLC_B{k}{l}(:);
                  opening_pos{k}(:,2) = MLC_A{k}{l}(:);
                  
                  %jaws:
                  opening_pos_x1{k}(1,1) = X1_JAW{k}{l};
                  opening_pos_x2{k}(1,1) = X2_JAW{k}{l};
                  opening_pos_y1{k}(1,1) = Y1_JAW{k}{l}+10;
                  opening_pos_y2{k}(1,1) = Y2_JAW{k}{l}+10;
               
                  X1_JAW_ROI(1).Position(1) = X1_JAW_ROI(1).Position(1)-(opening_pos_x1{k}(1,1));
                  X2_JAW_ROI(1).Position(1) = X2_JAW_ROI(1).Position(1)+(opening_pos_x2{k}(1,1));
                  Y1_JAW_ROI(1).Position(2) = Y1_JAW_ROI(1).Position(2)+(opening_pos_y1{k}(1,1));
                  Y2_JAW_ROI(1).Position(2) = Y2_JAW_ROI(1).Position(2)-(opening_pos_y2{k}(1,1));
        
               
                   for m = 1:length(MLC_A{k}{l})% MLC positions for one control point control point
        
                        
                            
                    %update current position to next position
                    A(m).Position(1) = A(m).Position(1)-(opening_pos{k}(m,1)); %leaf offset A
                    B(m).Position(1) = B(m).Position(1)+(opening_pos{k}(m,2)); %leaf offset A
                         
                     mask_A = ~createMask(A(m),fluence_grid);
                     mask_B = ~createMask(B(m),fluence_grid);
    
                     MLC_mask = MLC_mask.*(mask_A.*mask_B);

                   
                   end
             
                %MLC_mask_sum{k} = cellfun(@sum,MLC_mask,'UniformOutput',false)   
                %create tha jaw mask and update the fluence
                mask_X1 = ~createMask(X1_JAW_ROI(1),fluence_grid);
                mask_X2 = ~createMask(X2_JAW_ROI(1),fluence_grid);
                mask_Y1 = ~createMask(Y1_JAW_ROI(1),fluence_grid);
                mask_Y2 = ~createMask(Y2_JAW_ROI(1),fluence_grid);
                
                %create the masks
                JAW_mask{l} = mask_X1.*mask_X2.*mask_Y1.*mask_Y2;
                MLC_mask_sum{l} = MLC_mask;
                JAW_MLC_mask{l} = MLC_mask_sum{l}.*JAW_mask{l};
        
                %reset the MLC_mask
                MLC_mask = ones(400,400);
                %the accumulator
                %fluence_grid_total = (fluence_grid_total+(fluence_grid.*JAW_MLC_mask{l}))+MU_frac{k}(l);
                fluence_grid_total = (fluence_grid_total+(fluence_grid.*JAW_MLC_mask{l}.*MU_seg(l,k)));

                for i = 1:length(A)
                    A(i).Position(1) = MLC_A_shift(i);
                    B(i).Position(1) = MLC_B_shift(i);
                end
                
                %reset the JAW ROI
                X1_JAW_ROI(1).Position(1) = JAW_X1_shift(1);
                X2_JAW_ROI(1).Position(1) = JAW_X2_shift(1);
                Y1_JAW_ROI(1).Position(2) = JAW_Y1_shift(1);
                Y2_JAW_ROI(1).Position(2) = JAW_Y2_shift(1);
                
            
                %multiWaitbar( 'Processing Control Point(s)...', (l/length(MLC_A{k})), 'Color', 'g' );   
                if CollAngle(l,k)<180.1 %clockwise rotation
                %final_fluence_beam{k} = flipud(imrotate(fluence_grid_total,1*CollAngle(k,1),'bilinear','crop'));
                    final_fluence_beam_A{k} = (imrotate(fluence_grid_total,1*CollAngle(l,k),'bilinear','crop'));
                    %final_fluence_plan = imrotate(final_fluence_plan,-1*CollAngle(k,1),'bilinear','crop');
                else %counter clockwise rotation
                     %final_fluence_beam{k} = flipud(imrotate(fluence_grid_total,-1*CollAngle(k,1),'bilinear','crop'));
                     final_fluence_beam_A{k} = (imrotate(fluence_grid_total,-1*CollAngle(l,k),'bilinear','crop'));
                    %final_fluence_plan = imrotate(final_fluence_plan,1*CollAngle(k,1),'bilinear','crop');
                end
            end
       
        %add them to the final plan fluence and then reset the coll
        %final_fluence_beam{k} = flipud(fluence_grid_total);
        fluence_grid_total = zeros(400,400);
        final_fluence_plan_A = final_fluence_plan_A + final_fluence_beam_A{k};
        disp(strcat(OutFolder,fname(1,1:end-3),'.csv'))
        writematrix(final_fluence_beam_A{k},strcat(OutFolder,fname(1,:),'_',VarianSerial{k},'_',num2str(k),'_A','.csv'));
        %if (k > 1) && (k == size_subBeam)
        %        writematrix(final_fluence_plan_A,strcat(OutFolder,fname(1,1:end-3),'_',num2str(k),'_TOTAL','.csv'));
        %end
       
        
        fluence_grid_total = zeros(400,400);
        MLC_mask = ones(400,400);
        final_fluence_plan_E = zeros(400,400);

        
         if isempty(MLC_Ae{k})
                continue
         end   
         
            for l = 1:length(MLC_Ae{k})

                  if isempty(MLC_Be{k}{l})
                      continue
                  end
                  opening_pos{k}(:,1) = MLC_Be{k}{l}(:);
                  opening_pos{k}(:,2) = MLC_Ae{k}{l}(:);
                  
                  %jaws:
                  opening_pos_x1{k}(1,1) = X1_JAWe{k}{l};
                  opening_pos_x2{k}(1,1) = X2_JAWe{k}{l};
                  opening_pos_y1{k}(1,1) = Y1_JAWe{k}{l}+10;
                  opening_pos_y2{k}(1,1) = Y2_JAWe{k}{l}+10;
               
                  X1_JAW_ROI(1).Position(1) = X1_JAW_ROI(1).Position(1)-(opening_pos_x1{k}(1,1));
                  X2_JAW_ROI(1).Position(1) = X2_JAW_ROI(1).Position(1)+(opening_pos_x2{k}(1,1));
                  Y1_JAW_ROI(1).Position(2) = Y1_JAW_ROI(1).Position(2)+(opening_pos_y1{k}(1,1));
                  Y2_JAW_ROI(1).Position(2) = Y2_JAW_ROI(1).Position(2)-(opening_pos_y2{k}(1,1));
        
               
                   for m = 1:length(MLC_Ae{k}{l})% MLC positions for one control point control point
        
                        
                            
                    %update current position to next position
                    A(m).Position(1) = A(m).Position(1)-(opening_pos{k}(m,1)); %leaf offset A
                    B(m).Position(1) = B(m).Position(1)+(opening_pos{k}(m,2)); %leaf offset A
                         
                     mask_A = ~createMask(A(m),fluence_grid);
                     mask_B = ~createMask(B(m),fluence_grid);
    
                     MLC_mask = MLC_mask.*(mask_A.*mask_B);

                   
                   end
             
                %MLC_mask_sum{k} = cellfun(@sum,MLC_mask,'UniformOutput',false)   
                %create tha jaw mask and update the fluence
                mask_X1 = ~createMask(X1_JAW_ROI(1),fluence_grid);
                mask_X2 = ~createMask(X2_JAW_ROI(1),fluence_grid);
                mask_Y1 = ~createMask(Y1_JAW_ROI(1),fluence_grid);
                mask_Y2 = ~createMask(Y2_JAW_ROI(1),fluence_grid);
                
                %create the masks
                JAW_mask{l} = mask_X1.*mask_X2.*mask_Y1.*mask_Y2;
                MLC_mask_sum{l} = MLC_mask;
                JAW_MLC_mask{l} = MLC_mask_sum{l}.*JAW_mask{l};
        
                %reset the MLC_mask
                MLC_mask = ones(400,400);
                %the accumulator
                %fluence_grid_total = (fluence_grid_total+(fluence_grid.*JAW_MLC_mask{l}))+MU_frac{k}(l);
                fluence_grid_total = (fluence_grid_total+(fluence_grid.*JAW_MLC_mask{l}.*MU_sege(l,k)));

                for i = 1:length(A)
                    A(i).Position(1) = MLC_A_shift(i);
                    B(i).Position(1) = MLC_B_shift(i);
                end
                
                %reset the JAW ROI
                X1_JAW_ROI(1).Position(1) = JAW_X1_shift(1);
                X2_JAW_ROI(1).Position(1) = JAW_X2_shift(1);
                Y1_JAW_ROI(1).Position(2) = JAW_Y1_shift(1);
                Y2_JAW_ROI(1).Position(2) = JAW_Y2_shift(1);
                
            
                %multiWaitbar( 'Processing Control Point(s)...', (l/length(MLC_A{k})), 'Color', 'g' );   
                if CollAnglee(l,k)<180.1 %clockwise rotation
                %final_fluence_beam{k} = flipud(imrotate(fluence_grid_total,1*CollAngle(k,1),'bilinear','crop'));
                    final_fluence_beam_E{k} = (imrotate(fluence_grid_total,1*CollAnglee(l,k),'bilinear','crop'));
                    %final_fluence_plan = imrotate(final_fluence_plan,-1*CollAngle(k,1),'bilinear','crop');
                else %counter clockwise rotation
                     %final_fluence_beam{k} = flipud(imrotate(fluence_grid_total,-1*CollAngle(k,1),'bilinear','crop'));
                     final_fluence_beam_E{k} = (imrotate(fluence_grid_total,-1*CollAnglee(l,k),'bilinear','crop'));
                    %final_fluence_plan = imrotate(final_fluence_plan,1*CollAngle(k,1),'bilinear','crop');
                end
            end
       
        %add them to the final plan fluence and then reset the coll
        %final_fluence_beam{k} = flipud(fluence_grid_total);
            fluence_grid_total = zeros(400,400);
            final_fluence_plan_E = final_fluence_plan_E + final_fluence_beam_E{k};
            disp(strcat(OutFolder,fname(1,1:end-3),'.csv'))
            writematrix(final_fluence_beam_E{k},strcat(OutFolder,fname(1,:),'_',VarianSerial{k},'_',num2str(k),'_E','.csv'));

            %if (k > 1) && (k == size_subBeam)
            %    writematrix(final_fluence_plan_E,strcat(OutFolder,fname(1,1:end-3),'_',num2str(k),'_E_TOTAL','.csv'));
            %end
    

    if isempty(final_fluence_beam_E{k})
        continue
    end
    %perform the gamma analysis
    Image1 = final_fluence_beam_A{k};
    Image2 = final_fluence_beam_E{k};
    
    res_x = 1;
    res_y = 1;
    
    DTA_tol = GammaDTA
    radlim = DTA_tol * 1.5;
    rad = min(ceil(radlim/res_x),ceil(radlim/res_y));
    MaxVal = max(Image2(:));
    FE_thresh = GammaThresh;
    Dose_tol = GammaDD;
    %MaxVal = max(max(Image1));
    %[Im1_ydim, Im1_xdim] = size(Image1);
    Mask = zeros(size(Image2));
    crit_val = FE_thresh*MaxVal;
    
    %[Im1_ydim Im1_xdim] = size(Image1);
    
    % Use the Resampled image for this - EPID no spikes?
    Mask(Image2>crit_val) = 1;
    
    Dose_tol = Dose_tol*MaxVal; % GLOBAL - PERCENTAGE OF MAX DOSE
    
    % VECTORIZED CALCULATION STARTS HERE
    
    % GammaMapsub will carry the calculated gamma values for the truncated
    % images. GammaMap2 will be the Gamma values for the full image.
    GammaMapsub = NaN;
    GammaMap = zeros(size(Image1)).*Image1;  % casts to gpu if present
    
    % Find the threshold limits for truncation
    [validmask_y, validmask_x] = find(Mask);
    min_x = min(validmask_x)-rad;
    max_x = max(validmask_x)+rad;
    min_y = min(validmask_y)-rad;
    max_y = max(validmask_y)+rad;
    if min_x < 1
        min_x = 1;
    end
    if min_y < 1
        min_y = 1;
    end
    if max_x > size(Image1,2)
        max_x = size(Image1,2);
    end
    if max_y > size(Image1,1)
        max_y = size(Image1,1);
    end
    % Truncate the images to avoid needless calculations
    Im1 = Image1(min_y:max_y,min_x:max_x);
    Im2 = Image2(min_y:max_y,min_x:max_x);
    % Shift the image by varying amounts. Determine the minimum gamma value
    % for all shifts
    l = 1;
    for i=-rad:rad
        for j=-rad:rad
            % circshift function wraps elements from top to bottom as necessary
            % The entire image is shifted at once
            Im2_shift = circshift(Im2,[i j]);
            dist = sqrt((res_y*i)^2 + (res_x*j)^2);
            DoseDiff = Im2_shift - Im1;
            % Compute the gamma map for this particular shift value
            Gamma_temp = sqrt((dist./DTA_tol).^2 + (DoseDiff./Dose_tol).^2);
            %Gamma_temp = sqrt((dist./DTA_tol).^2 + (DoseDiff).^2);
            % Accumulate the map of the minimum values of gamma at each point
            GammaMapsub = min(GammaMapsub,Gamma_temp);
        end
    l = l+1;
    
    end
    
    
    % Put the truncated gamma map back into its proper location within the full
    % gamma map
    GammaMap(min_y:max_y,min_x:max_x) = GammaMapsub;
    % Remove any edge effects from the circular shifting by multiplying by the
    % mask values. This will negate any calculated gamma values around the
    % edges of the distribution where this effect would arise
    GammaMapNoMask = GammaMap;
    GammaMap = GammaMap .* Mask;
    % Ensure that NaN values outside the mask do not affect the calculation
    GammaMap(~Mask) = 0.0;
    
    % Compute statistics
    numWithinField = nnz(Mask);
    numpass = nnz(GammaMap<1 & Mask)./numWithinField;
    avg = sum(GammaMap(:))./numWithinField;
    
    %unpad the array
    pp_init{k,1} = numpass*100; %%
    
    if pp_init{k,1} > (PassingCriteria*100)
        Result{k,1} = 'PASS';
    else
        Result{k,1} = 'FAIL';
    end
    pause(5)

    %set the fielname formats
    writematrix(GammaMap,strcat(OutFolder,fname(1,:),'_',VarianSerial{k},'_',num2str(k),'_G','.csv'));
    BEAM_ON{k,1} = BeamOnTime(k,1);
    fname_A{k,1} = strcat(OutFolder,fname(1,:),'_',VarianSerial{k},'_',num2str(k),'_A','.csv');
    fname_E{k,1} = strcat(OutFolder,fname(1,:),'_',VarianSerial{k},'_',num2str(k),'_E','.csv');
    fname_G{k,1} = strcat(OutFolder,fname(1,:),'_',VarianSerial{k},'_',num2str(k),'_G','.csv');
    % = fname_out;
     %%CALCULATE ALL OUR ERRORS
    RMS_MLC{k,1} = max(sqrt(mean([snapShotFinal.MLC_A]-[snapShotFinal.MLC_E]).^2));
    MAX_X1{k,1} = max([snapShotFinal.X1_A] - [snapShotFinal.X1_E]);
    MAX_X2{k,1} = max([snapShotFinal.X2_A] - [snapShotFinal.X2_E]);
    MAX_Y1{k,1} = max([snapShotFinal.Y1_A] - [snapShotFinal.Y1_E]);
    MAX_Y2{k,1} = max([snapShotFinal.Y2_A] - [snapShotFinal.Y2_E]);
    MAX_LAT{k,1} = max([snapShotFinal.couchLatA] - [snapShotFinal.couchLatE]);
    MAX_LONG{k,1} = max([snapShotFinal.couchLngA] - [snapShotFinal.couchLngE]);
    MAX_VERT{k,1} = max([snapShotFinal.couchVrtA] - [snapShotFinal.couchVrtE]);
    MAX_COLL{k,1} = max([snapShotFinal.colRotationA] - [snapShotFinal.colRotationE]);
    MAX_GANTRY{k,1} = max([snapShotFinal.gantryRotationA] - [snapShotFinal.gantryRotationE]);
    ID{k,1} = extractBefore(fname_out,'_');
    Plan{k,1} = extractAfter(fname_out,'_');
    Beam{k,1} = extractBetween(Plan(k,1),'_','_');
    Beam{k,1} = Beam{k,1}(end)
    Beam{k,1} = strcat(Beam{k,1},'_',num2str(k));
    %set the date format
    Y = fname_out(1,end-13:end-13+3); 
    M = fname_out(1,end-9:end-10+2);
    D = fname_out(1,end-7:end-7+1);
    H = fname_out(1,end-5:end-5+1);
    MI = fname_out(1,end-3:end-3+1);
    S = fname_out(1,end-1:end-1+1);
    LogDate{k,1} = [Y,'-',M,'-',D,' ', H,':',MI,':',S];

    
end

end