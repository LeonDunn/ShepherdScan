function eFluenceQA_Agility(fileName, OutFolder)
%get the number of films and
if iscell(fileName)
    num_films = length(fileName);
    for q=1:num_films   
        dcm_data{q,1} = dicomread(fileName);
        dcm_info{q,1} = dicominfo(fileName);
        PtName{q,1} = dcm_info{1}.PatientName;
        PtID{q,1} = dcm_info{1}.PatientID;
        PtDOB{q,1} = dcm_info{1}.PatientBirthDate;
        PtSEX{q,1} = dcm_info{1}.PatientSex;
        Beams{q,1} = struct2cell(dcm_info{1}.BeamSequence);
        TxMACHINE{q,1} = Beams{1}{1}.TreatmentMachineName;
        Plan{q,1} = dcm_info{1}.RTPlanLabel;

    end
else
    num_films = 1;
    for q=1:num_films   
        dcm_data{q,1} = dicomread(fileName);
        dcm_info{q,1} = dicominfo(fileName);
        PtName{q,1} = dcm_info{1}.PatientName;
        PtID{q,1} = dcm_info{1}.PatientID;
        PtDOB{q,1} = dcm_info{1}.PatientBirthDate;
        PtSEX{q,1} = dcm_info{1}.PatientSex;
        Beams{q,1} = struct2cell(dcm_info{1}.BeamSequence);
        TxMACHINE{q,1} = Beams{1}{1}.TreatmentMachineName;
        Plan{q,1} = dcm_info{1}.RTPlanLabel;

    end%
end



%for each beam in the baems, get the gantry angle, MU total etc.
    
    num_beams = size(Beams{1,1},1);
    j = 1;
    for i = 1:num_beams
    
        %get the beam type and only move on if its' dynamic
        
        if string(Beams{1,1}{i,1}.BeamType) == "STATIC"
            continue
        end
        %structs of structs
        CPS{j} = struct2cell(Beams{1,1}{i,1}.ControlPointSequence);
        BeamDesc{j} = Beams{1,1}{i,1}.BeamDescription
        MU_Total{j} = struct2cell(dcm_info{1}.FractionGroupSequence.Item_1.ReferencedBeamSequence);
        j = j+1;
    end
    
    %think this only works for jaw tracking...
    for i = 1:size(CPS,2)
        CollAngle(i,1) = CPS{1,i}{1,1}.BeamLimitingDeviceAngle;
        
        for j = 1:size(CPS{1,i},1)
            size(CPS{1,i},1);
            MU_frac{i}(j,1) = CPS{1,i}{j,1}.CumulativeMetersetWeight;
            fields = fieldnames(CPS{1,i}{j,1}.BeamLimitingDevicePositionSequence);
            if size(fields,1)>1 %must contain jaws
            
%                 X1_JAW{i}{j} = CPS{1,i}{j,1}.BeamLimitingDevicePositionSequence.Item_1.LeafJawPositions(1);
%                 X2_JAW{i}{j} = CPS{1,i}{j,1}.BeamLimitingDevicePositionSequence.Item_1.LeafJawPositions(2);
                Y1_JAW{i}{j} = CPS{1,i}{j,1}.BeamLimitingDevicePositionSequence.Item_1.LeafJawPositions(1);
                Y2_JAW{i}{j} = CPS{1,i}{j,1}.BeamLimitingDevicePositionSequence.Item_1.LeafJawPositions(2);
                MLC = CPS{1,i}{j,1}.BeamLimitingDevicePositionSequence.Item_2.LeafJawPositions;
                MLC_A{i}{j} = flipud(MLC(1:80,1));
                MLC_B{i}{j} = flipud(MLC(81:160,1));
                MU{i}{j} = CPS{1,i}{j,1}.ReferencedDoseReferenceSequence.Item_1.CumulativeDoseReferenceCoefficient;
                MU_Field{i}{j} = MU_Total{1,1}{i,1}.BeamMeterset;
                MU_seg(j,i) = MU{i}{j}*MU_Field{i}{j};
                
            else
%               X1_JAW{i}{j} = CPS{1,i}{1,1}.BeamLimitingDevicePositionSequence.Item_1.LeafJawPositions(1);
%               X2_JAW{i}{j} = CPS{1,i}{1,1}.BeamLimitingDevicePositionSequence.Item_1.LeafJawPositions(2);
                Y1_JAW{i}{j} = CPS{1,i}{1,1}.BeamLimitingDevicePositionSequence.Item_1.LeafJawPositions(1);
                Y2_JAW{i}{j} = CPS{1,i}{1,1}.BeamLimitingDevicePositionSequence.Item_2.LeafJawPositions(2);
                
                MLC = CPS{1,i}{j,1}.BeamLimitingDevicePositionSequence.Item_2.LeafJawPositions;
                MLC_A{i}{j} = flipud(MLC(1:80,1));
                MLC_B{i}{j} = flipud(MLC(81:160,1));
                MU{i}{j} = CPS{1,i}{j,1}.ReferencedDoseReferenceSequence.Item_1.CumulativeDoseReferenceCoefficient;
                MU_Field{i}{j} = MU_Total{1,1}{i,1}.BeamMeterset;
                MU_seg(j,i) = MU{i}{j}*MU_Field{i}{j};
                
            end
            
        end
    end
            
    %correct the coll angles that are counter clockwise
    MU_seg = [repmat(0,[1,size(CPS,2)]); diff(MU_seg)];
    MU_seg(MU_seg<0) = 0;
    %start at zero and accumlate fluence
    fluence_grid = ones(400,400);
    fluence_grid_total = zeros(400,400);
    MLC_mask = ones(400,400);
    final_fluence_plan = zeros(400,400);

    grid on
    axis on
    %n_large_leaves = 10;
    n_small_leaves = 80;
    leaf_width_large = 5; %mm
    leaf_width_small = 5; %mm
    large_leaf_offset = leaf_width_large;
    leaf_length = 200;
    num_leaves = 80;
    
        
    
    %no graphics
    j = 1;
    %create the rectangle MLC leaves
    for i = 0:5:395
        A(j) = images.roi.Rectangle('Position',[0,i,200,leaf_width_small],'Color','r','InteractionsAllowed','none');
        %A(j).FaceAlpha = 0.1;
        j=j+1;
    end
        
    j = 1;
    %create the rectangle MLC leaves
    for i = 0:5:395
        B(j) = images.roi.Rectangle('Position',[200,i,200,leaf_width_large],'Color','b','InteractionsAllowed','none');
        B(j).FaceAlpha = 0.1;
        j=j+1;
    end
%     %create the rectangle MLC leaves
   
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

    
%     multiWaitbar( 'Processing Beam(s)...', 0, 'Color', 'g' );
%     multiWaitbar( 'Processing Control Point(s)...', 0, 'Color', 'g' );
%     multiWaitbar( 'Processing MLC Pair...', 0, 'Color', 'g' );
    
    
    
    for k = 1:length(MLC_A)
           
    
        for l = 1:length(MLC_A{k})
              
              opening_pos{k}(:,1) = MLC_A{k}{l}(:);
              opening_pos{k}(:,2) = MLC_B{k}{l}(:);
%               figure(101)
%               hold on
%               plot(opening_pos{k}(:,2))  

              %these are the max MLC positions used to create the final
              %mask for the most extended leaf
              A_max = min(opening_pos{k}(:,1));
              B_max = max(opening_pos{k}(:,2));
              %basically simlulate the x jaws ends by making anthing before the most extended a leaf zero and anything after the most extended b leaf zero  
%               fluence_grid(:,1:A_max) = 0;
%               fluence_grid(:,B_max:end) = 0;
              %jaws:
              opening_pos_x1{k}(1,1) = A_max-5;
              opening_pos_x2{k}(1,1) = B_max+5;
              opening_pos_y1{k}(1,1) = Y1_JAW{k}{l};
              opening_pos_y2{k}(1,1) = -1*Y2_JAW{k}{l};
              
              %this is setting y jaws backside to 0
              fluence_grid(1:floor(opening_pos_y2{k}(1,1)),:) = 0;
              fluence_grid(200-floor(opening_pos_y1{k}(1,1)):end,:) = 0;

              X1_JAW_ROI(1).Position(1) = X1_JAW_ROI(1).Position(1)+(opening_pos_x1{k}(1,1));
              X2_JAW_ROI(1).Position(1) = X2_JAW_ROI(1).Position(1)+(opening_pos_x2{k}(1,1));
              Y1_JAW_ROI(1).Position(2) = Y1_JAW_ROI(1).Position(2)-(opening_pos_y2{k}(1,1));
              Y2_JAW_ROI(1).Position(2) = Y2_JAW_ROI(1).Position(2)+(opening_pos_y1{k}(1,1));
              %Y1_JAW_ROI(1).Position(2) = Y1_JAW_ROI(1).Position(2)-(opening_pos_y1{k}(1,1));
              %Y2_JAW_ROI(1).Position(2) = Y2_JAW_ROI(1).Position(2)-(opening_pos_y2{k}(1,1)); 
           
               for m = 1:length(MLC_A{k}{l})% MLC positions for one control point control point
    
                      
                    %update current position to next position
                    A(m).Position(1) = A(m).Position(1)+(opening_pos{k}(m,1)); %leaf offset A
                    B(m).Position(1) = B(m).Position(1)+(opening_pos{k}(m,2)); %leaf offset A
                    

                     mask_A = ~createMask(A(m),fluence_grid);
                     mask_B = ~createMask(B(m),fluence_grid);
    
                     MLC_mask = MLC_mask.*(mask_A.*mask_B);
                     
                    
%                multiWaitbar( 'Processing MLC Pair...', (m/length(MLC_A{k}{l})), 'Color', 'g');
               %drawnow
               
               end
            %drawnow
            %return
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
            %fluence_grid_total = (fluence_grid_total+(fluence_grid.*JAW_MLC_mask{l}));
            %imshow(fluence_grid_total)
            %fluence_grid_total(fluence_grid_total~=0) = fluence_grid_total(fluence_grid_total~=0)+MU{k}{l};
            %csvwrite(strcat('Beam_',num2str(k),'_',num2str(l),'.csv'),fluence_grid_total)
            %imshow(fluence_grid_total)
            %imshow(double(cell2mat(MLC_mask_sum{k})))
    %        set the original positions to reset to
    %         if l ==2
    %                 return
    %         end
            for i = 1:length(A)
                A(i).Position(1) = MLC_A_shift(i);
                B(i).Position(1) = MLC_B_shift(i);
            end
            
            %reset the JAW ROI
            X1_JAW_ROI(1).Position(1) = JAW_X1_shift(1);
            X2_JAW_ROI(1).Position(1) = JAW_X2_shift(1);
            Y1_JAW_ROI(1).Position(2) = JAW_Y1_shift(1);
            Y2_JAW_ROI(1).Position(2) = JAW_Y2_shift(1);
        %reset the fluence grid back to one    
            fluence_grid = ones(400,400);
%         multiWaitbar( 'Processing Control Point(s)...', (l/length(MLC_A{k})), 'Color', 'g' );   
        end
    
    %rotate the collimator back to 0
    %fluence_grid = imrotate(fluence_grid,-1*CollAngle(k,1));
    %fluence_grid_total = imrotate(fluence_grid_total,-1*CollAngle(k,1));
    %MLC_mask = imrotate(MLC_mask,-1*CollAngle(k,1));
    %final_fluence_plan = imrotate(final_fluence_plan,-1*CollAngle(k,1)); 
    %this is the sum for the beam
    if CollAngle(k,1)<180.1 %clockwise rotation
        %final_fluence_beam{k} = flipud(imrotate(fluence_grid_total,1*CollAngle(k,1),'bilinear','crop'));
        final_fluence_beam{k} = (imrotate(fluence_grid_total,1*CollAngle(k,1),'bilinear','crop'));
        %final_fluence_plan = imrotate(final_fluence_plan,-1*CollAngle(k,1),'bilinear','crop');
    else %counter clockwise rotation
        %final_fluence_beam{k} = flipud(imrotate(fluence_grid_total,-1*CollAngle(k,1),'bilinear','crop'));
        final_fluence_beam{k} = (imrotate(fluence_grid_total,-1*CollAngle(k,1),'bilinear','crop'));
        %final_fluence_plan = imrotate(final_fluence_plan,1*CollAngle(k,1),'bilinear','crop');
    end
    %add them to the final plan fluence and then reset the coll
    final_fluence_beam{k} = (fluence_grid_total);
    fluence_grid_total = zeros(400,400);
    final_fluence_plan = final_fluence_plan + final_fluence_beam{k};
    figure(k)
    imagesc(final_fluence_beam{k})
    hold on
    x = [0 400];
    y = [200 200];
    pl = line(x,y,'Color','red','LineStyle','--','LineWidth',1.5);
    x = [200 200];
    y = [0 400];
    p2 = line(x,y,'Color','red','LineStyle','--','LineWidth',1.5);
    %set(gca,'YDir','reverse');
   
    writematrix(final_fluence_beam{k},strcat(OutFolder,PtID{1},'_',Plan{1},'_Beam_',num2str(k),'_',TxMACHINE{1},'.csv'));
    
%     csvwrite(strcat('Beam_',num2str(k),'_',num2str(CollAngle(k)),'.csv'),final_fluence_beam{k})
    % if CollAngle(k,1)<180.1 %clockwise rotation
    %     final_fluence_beam{k} = imrotate(fluence_grid_total,1*CollAngle(k,1),'bilinear','crop');
    %     %final_fluence_plan = imrotate(final_fluence_plan,1*CollAngle(k,1),'bilinear','crop');
    % else %counter clockwise rotation
    %     final_fluence_beam{k} = imrotate(fluence_grid_total,-1*CollAngle(k,1),'bilinear','crop');
    %     %final_fluence_plan = imrotate(final_fluence_plan,-1*CollAngle(k,1),'bilinear','crop');
    % end
    %reset the beam fluence
%     multiWaitbar( 'Processing Beam(s)...', (k/length(MLC_A)), 'Color', 'g' );
    
    end
    figure(k+1)
    imagesc(final_fluence_plan)
    writematrix(final_fluence_plan,strcat(OutFolder,PtID{1},'_',Plan{1},'_',TxMACHINE{1},'_TOTAL','.csv'))
end   
