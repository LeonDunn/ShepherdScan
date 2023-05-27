
%  close all
%  clear all
  filename =  '134254_N30 L Neck_1.1T_20230518100614.bin'
  OutFolder = 'C:\Temp\Output\'

    [folder, fname] = fileparts(filename);
    
    
    fid=fopen(filename);
    %
    
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
       fseek(fid,1024+(i*560),'bof'); %uncommment this line for version 2.0
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
    
    
    
    for i=1:numberOfSnapShots
        % actual_(1,j,a1),expected_(1,j,a1),actual_(1,j,a2),expected_(1,j,a2),....
        % see figure on page 11 of Specification
    
          %temp5{i}= single(fread(fid,2*totalNumberOfSamples,'single'));
    
          temp5{i} = single(fread(fid,2*totalNumberOfSamples,'*single'));
          %temp5{i} = double(temp5{i});
          BH = int32(temp5{i}(28));
          
          %if the beam is held, we just want to skip recording the data and
          %move on
          if (BH > 1) 
              %disp('beam holding')
              temp5{i} = [];
              continue
          end
          
         
    % need to convert angles (colli, gantry, couch) from intrinsic Varian scale (+1°) to IEC 61217 scale (+179°)
     
    % snapShot(i).colRotationE=single(temp5{i}(1));
          snapShot(i).colRotationE = Varian2IEC((temp5{i}(1)));
          
          snapShot(i).colRotationA = Varian2IEC((temp5{i}(2)));
        
          snapShot(i).gantryRotationE = Varian2IEC((temp5{i}(3)));
          
          snapShot(i).gantryRotationA = Varian2IEC((temp5{i}(4)));
          
          snapShot(i).Y1_E = (temp5{i}(5));
          
          snapShot(i).Y1_A = (temp5{i}(6));
          
          snapShot(i).Y2_E = (temp5{i}(7));
          
          snapShot(i).Y2_A = (temp5{i}(8));
          
          snapShot(i).X1_E = (temp5{i}(9));
          
          snapShot(i).X1_A = (temp5{i}(10));
          
          snapShot(i).X2_E = (temp5{i}(11));
          
          snapShot(i).X2_A = (temp5{i}(12));
          
          snapShot(i).couchVrtE = (temp5{i}(13));
          
          snapShot(i).couchVrtA = (temp5{i}(14));
          
          snapShot(i).couchLngE = (temp5{i}(15));
          
          snapShot(i).couchLngA = (temp5{i}(16));
          
          snapShot(i).couchLatE = (temp5{i}(17));
          
          snapShot(i).couchLatA = (temp5{i}(18));
          
          snapShot(i).couchRotationE = Varian2IEC((temp5{i}(19)));
          
          snapShot(i).couchRotationA = Varian2IEC((temp5{i}(20)));
    
          snapShot(i).MU_E = (temp5{i}(25));
          
          snapShot(i).MU_A = (temp5{i}(26));
          
          snapShot(i).beamHoldE = (temp5{i}(27));
          
          snapShot(i).beamHoldA = (temp5{i}(28));
          
          snapShot(i).controlPointE = (temp5{i}(29));
          
          snapShot(i).controlPointA = (temp5{i}(30));
          
          snapShot(i).CarrA_E = (temp5{i}(31));
          
          snapShot(i).CarrA_A = (temp5{i}(32));
          
          snapShot(i).CarrB_E = (temp5{i}(33));
          
          snapShot(i).CarrB_A = (temp5{i}(34));
          
          %last but not least comes the MLC
          
          snapShot(i).MLC_E = (temp5{i}(35:2:end-1));
          snapShot(i).MLC_A = (temp5{i}(36:2:end));
    
          
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
    i = 1;
      
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
                        disp('new beam')
                        
                            disp(j)
                            i = i+1;
                    end
                    
                                  
                    CollAngle(j,i) = snapShotFinal(j).colRotationA;
                    
                    X1_JAW{i}{j} = snapShotFinal(j).X1_A*10;
                    X2_JAW{i}{j} = snapShotFinal(j).X2_A*10;
                    Y1_JAW{i}{j} = snapShotFinal(j).Y1_A*10;
                    Y2_JAW{i}{j} = snapShotFinal(j).Y2_A*10;
                   
                    MLC_A{i}{j} = flipud(snapShotFinal(j).MLC_A(1:60,1)*10);
                    MLC_B{i}{j} = flipud(snapShotFinal(j).MLC_A(61:end,1)*10);
                    MU_seg(j,i) = snapShotFinal(j).MU_A;
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
        end
    
    end
    
   
    
    %correct the coll angles that are counter clockwise
    MU_seg = [repmat(0,[1,size(MU_seg,2)]); diff(MU_seg)];
    MU_seg(MU_seg<0) = 0;
    MU_seg(MU_seg>100) = 0;
    
    fluence_grid = ones(400,400);
 
        
    grid on
    axis on
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
        final_fluence_plan = zeros(400,400);
    
        for k = 1:length(MLC_A)

         
            for l = 1:length(MLC_A{k})
    
             %for l = 3450:length(MLC_A{k})
                  if isempty(MLC_B{k}{l})
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
                    final_fluence_beam{k} = (imrotate(fluence_grid_total,1*CollAngle(l,k),'bilinear','crop'));
                    %final_fluence_plan = imrotate(final_fluence_plan,-1*CollAngle(k,1),'bilinear','crop');
                else %counter clockwise rotation
                     %final_fluence_beam{k} = flipud(imrotate(fluence_grid_total,-1*CollAngle(k,1),'bilinear','crop'));
                     final_fluence_beam{k} = (imrotate(fluence_grid_total,-1*CollAngle(l,k),'bilinear','crop'));
                    %final_fluence_plan = imrotate(final_fluence_plan,1*CollAngle(k,1),'bilinear','crop');
                end
            end
       
        %add them to the final plan fluence and then reset the coll
        %final_fluence_beam{k} = flipud(fluence_grid_total);
        fluence_grid_total = zeros(400,400);
        final_fluence_plan = final_fluence_plan + final_fluence_beam{k};
        disp(strcat(OutFolder,fname(1,1:end-3),'.csv'))
        writematrix(final_fluence_beam{k},strcat(OutFolder,fname(1,1:end-3),'_',num2str(k),'.csv'));
   
        end
              

