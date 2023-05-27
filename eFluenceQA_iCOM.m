function [Name, ID, Field, Machine, Rx] = eFluenceQA_iCOM(fileNameiCOM, OutFolder)
%get the number of films and
icom_chk = (strfind(fileNameiCOM(1,:),'.xz'));
                
if ~isempty(icom_chk) > 0 %found a iCOM
    %we need to unzip it
    UnzipiCOM(fileNameiCOM)
    %now the filename is simply the filename minus the .xz
    %bit
    fileNameiCOM = fileNameiCOM(1,1:end-3);

end

[fPath,fName,ext] = fileparts(fileNameiCOM)

disp('im runnnin')
    %get the number of films and
if iscell(fileNameiCOM)
    num_films = length(fileNameiCOM);
else
    num_films = 1;
end

        IDExp = '\x00\LO\x00\P\x06\x00\x00\x00[a-zA-Z_,0-9 ]*';
        FieldExp = '\x03\x10\LO\x00\P\x0E\x00\x00\x00[.a-zA-Z_,0-9 \/\\\-]*';
        FieldExp2 = '^\x00\x00\x00[.a-zA-Z_,0-9 \/\\\-]*';    
        RxExp = '\x01\x10\LO\x00\P(\x14|\x0E|\x0f|\x08)\x00\x00\x00[.a-zA-Z_,0-9 \/\\\-]*';
        RxExp2 = '^\x00\x00\x00[.a-zA-Z_,0-9 \/\\\-]*\x01';
        BeamExp = '\x03\x10\LO\x00\P(\x0E|\x0B)\x00\x00\x00[.a-zA-Z_,0-9 \/\\\-]*';
        BeamExp2 = '^\x00\x00\x00[.a-zA-Z_,0-9 \/\\\-]*';
        
        for i=1:num_films %loop through files
    
                    %check if there is more than one files then read in the
                    %PD data in terms of the metadata dicominfo and the
                    %actual data. We need the info as it contains the
                    %conversion information for dose and pixel spacing
               
        iCOM_data{i,1} = importdata(fileNameiCOM);
        iCOM_info{i,1} = 'nothing';

        
        

        IDMatches = regexpi(iCOM_data{i},IDExp,'match', 'once');
        IDMatches = IDMatches(~cellfun(@isempty, IDMatches));
        IDMatches = regexprep(IDMatches,'[^\/.a-zA-Z0-9]',' ');   
        ID{i} = IDMatches

        FieldMatches = regexpi(iCOM_data{i},FieldExp,'match', 'once');
        FieldMatches = FieldMatches(~cellfun(@isempty, FieldMatches));
        FieldMatches = regexprep(FieldMatches,'[^\-\/.a-z-A-Z0-9]','');
        Field{i} = FieldMatches; %remove the LOP

        if isempty(Field{i})
            %try the second one
            FieldMatches = regexpi(iCOM_data{i},FieldExp2,'match', 'once');
            FieldMatches = FieldMatches(~cellfun(@isempty, FieldMatches));
            FieldMatches = regexprep(FieldMatches,'[^\-\/.a-z-A-Z0-9]','');
            Field{i} = FieldMatches; %remove the LOP
        end

        RxMatches = regexpi(iCOM_data{1},RxExp,'match', 'once');
        RxMatches = RxMatches(~cellfun(@isempty, RxMatches));
        RxMatches = regexprep(RxMatches,'[^ \,\-\.a-z-A-Z0-9]','');
        RxMatches = regexprep(RxMatches,'LOP','');
        Rx{i} = RxMatches;

        if isempty(Rx{i})
            %try the second one
            RxMatches = regexpi(iCOM_data{i},RxExp2,'match', 'once');
            RxMatches = RxMatches(~cellfun(@isempty, RxMatches));
            RxMatches = regexprep(RxMatches,'[^\-\/.a-z-A-Z0-9]','');
            
            Rx{i} = RxMatches; %remove the LOP
        end

        BeamMatches = regexpi(iCOM_data{1},BeamExp,'match', 'once');
        BeamMatches = BeamMatches(~cellfun(@isempty, BeamMatches));
        BeamMatches = regexprep(BeamMatches,'[^ \,\-\.a-z-A-Z0-9]','');
        BeamMatches = regexprep(BeamMatches,'LOP','');
        Beam{i} = BeamMatches

        if isempty(Beam{i})
            %try the second one
            BeamMatches = regexpi(iCOM_data{i},BeamExp2,'match', 'once');
            BeamMatches = BeamMatches(~cellfun(@isempty, BeamMatches));
            BeamMatches = regexprep(BeamMatches,'[^ \,\-\.a-z-A-Z0-9]','');
            Beam{i} = BeamMatches; %remove the LOP
        end

        
     
        %unique(Field)
        end
%begion processing
        j = 1;

 %set up the regexp
        %regular expressions for the JAMs and MLC - delivered
        dMLC_AExp = '\x1c\x01\DS\x00\R\x05\x00\x00\x00[-]?([1-9]{1}[0-9]{0,}(\.[0-9]{0,2})?|0(\.[0-9]{0,2})?|\.[0-9]{1,2})';
        dMLC_BExp = '\x1c\x01\DS\x00\R(\x04|\x06)\x00\x00\x00[-]?([1-9]{1}[0-9]{0,}(\.[0-9]{0,2})?|0(\.[0-9]{0,2})?|\.[0-9]{1,2})';
        dMLC_Exp =  '\x1c\x01\DS\x00\R(\x05|\x04|\x06)\x00\x00\x00[-]?([1-9]{1}[0-9]{0,}(\.[0-9]{0,2})?|0(\.[0-9]{0,2})?|\.[0-9]{1,2})';
        %regular expressions for the JAW and MLC this needs to be rearranged to
        %create the planned fluence by 
        %pMLC_XExp = '\x1c\x01\DS\x00\S\x05\x00\x00\x00[\.a-zA-Z_0-9\-]*'
        pMLC_Exp = '\x1c\x01\DS\x00\S(\x05|\x04|\x06)\x00\x00\x00[-]?([1-9]{1}[0-9]{0,}(\.[0-9]{0,2})?|0(\.[0-9]{0,2})?|\.[0-9]{1,2})';
        
        %get the field ID to sort the MLC by and other parameters
        NameExp = '\x00\PN\x00\P\x12\x00\x00\x00[.a-zA-Z_,0-9 \/\\\-]*';
       
        
        
        MachineExp = '\xb2\x00\SH\x00P\x04\x00\x00\x00[.a-zA-Z_,0-9 \/\\\-]*';
        RadiationExp = '0\xC6\x00\w.\x00R\x06\x00\x00\x00[.a-zA-Z_,0-9 \/\\\-]*';
        EnergyExp = '0\x14\x01\w.\x00\R.\x00\x00\x00[.a-zA-Z_,0-9 \/\\\-]*';
        SegmentExp = '\x07\x10\w.\x00\R\x01\x00\x00\x00\d';
        TotalMUExp = '\t\x10\DS\x00\R(\x05|\x04|\x03)\x00\x00\x00[-]?([1-9]{1}[0-9]{0,}(\.[0-9]{0,2})?|0(\.[0-9]{0,2})?|\.[0-9]{1,2})';
        DeliveryMUExp = '\x08\d+\x00\DS\x00\R(\x05|\x04|\x03)\x00\x00\x00[-]?([1-9]{1}[0-9]{0,}(\.[0-9]{0,2})?|0(\.[0-9]{0,2})?|\.[0-9]{1,2})'; %fixed
        %BackupDeliveryMUExp = '\x08\\x00\DS\x00\R\x03\x00\x00\x00[.a-zA-Z_,0-9 \/\\\-]*'
        BeamTimerExp = '8\x00\SH\x00\R(\x05|\x04|\x03)\x00\x00\x00[.a-zA-Z_,0-9 \/\\\-]*';
        SegmentMUExp = '\x0b\x00\DS\x00\R(\x05|\x04|\x03)\x00\x00\x00[.a-zA-Z_,0-9 \/\\\-]*';
        GantryExp = '\x1e\x01\DS\x00\R(\x03|\x05|\x04|\x06)\x00\x00\x00[-]?([1-9]{1}[0-9]{0,}(\.[0-9]{0,2})?|0(\.[0-9]{0,2})?|\.[0-9]{1,2})' ;
        CollimatorExp = ' \x01\DS\x00\R(\x03|\x05|\x04|\x06)\x00\x00\x00[-]?([1-9]{1}[0-9]{0,}(\.[0-9]{0,2})?|0(\.[0-9]{0,2})?|\.[0-9]{1,2})';
        BeamDescriptionExp = '\x00\SH\x00\R\x00\x00\x00[a-zA-Z_0-9\-]*';
        
        
        EnergyMatches = regexpi(iCOM_data{1},EnergyExp,'match', 'once');
        EnergyMatches = EnergyMatches(~cellfun(@isempty, EnergyMatches));
        EnergyMatches = regexprep(EnergyMatches,'[^\-\.a-z-A-Z0-9]','');
        Energy{j} = EnergyMatches;
        
        NameMatches = regexpi(iCOM_data{1},NameExp,'match', 'once');
        NameMatches = NameMatches(~cellfun(@isempty, NameMatches));
        NameMatches = regexprep(NameMatches,'[^\-\.a-z-A-Z0-9\,]','');
        Name{j} = NameMatches;
        
        
        
        MachineMatches = regexpi(iCOM_data{1},MachineExp,'match', 'once');
        MachineMatches = MachineMatches(~cellfun(@isempty, MachineMatches));
        MachineMatches = regexprep(MachineMatches,'[^\-\.a-z-A-Z0-9]','');
        Machine{j} = MachineMatches;

        RadiationMatches = regexpi(iCOM_data{1},RadiationExp,'match', 'once');
        RadiationMatches = RadiationMatches(~cellfun(@isempty, RadiationMatches));
        RadiationMatches = regexprep(RadiationMatches,'[^\-\.a-z-A-Z0-9]','');
        Radiation{j} = RadiationMatches;

        SegmentMatches = regexpi(iCOM_data{1},SegmentExp,'match','once');
        SegmentMatches = SegmentMatches(~cellfun(@isempty, SegmentMatches));
        SegmentMatches = regexprep(SegmentMatches,'[^\-\.0-9]','');
        SegmentMatches{j} = str2double(SegmentMatches);

        TotalMUMatches = regexpi(iCOM_data{1},TotalMUExp,'match','once');
        TotalMUMatches = TotalMUMatches(~cellfun(@isempty, TotalMUMatches));
        TotalMUMatches = regexprep(TotalMUMatches,'[^\-\.0-9]','');
        TotalMU{j} = str2double(TotalMUMatches);
        %calculate the difference in the total MU - use this as a skip feature
        diffMU{j} = diff(TotalMU{j});
        diffMU{j}(end+1) = 0;
        
       
        
        DeliveryMUMatches = regexpi(iCOM_data{1},DeliveryMUExp,'match','once');
        DeliveryMUMatches = DeliveryMUMatches(~cellfun(@isempty, DeliveryMUMatches));
        DeliveryMUMatches = regexprep(DeliveryMUMatches,'02','');
        DeliveryMUMatches = regexprep(DeliveryMUMatches,'[^\-\.0-9]','');
        DeliveryMU{j} = str2double(DeliveryMUMatches);

        
         %now we need to remove irrelevant control points
        %take the delivery MU and determine whether it's changing 
        diffDelMU{j} = [0; diff(DeliveryMU{j})];
        
        for k = 1: length(diffDelMU)
            ind_cp_mu{k} = find(diffDelMU{k}<0);
            diffDelMU{k}(ind_cp_mu{k}) = 0;
            ind_cp_mu{k} = find(diffDelMU{k}>100);
            diffDelMU{k}(ind_cp_mu{k}) = 0;
        end
        

        
        SegmentMUMatches = regexpi(iCOM_data{1},SegmentMUExp,'match','once');
        SegmentMUMatches = SegmentMUMatches(~cellfun(@isempty, SegmentMUMatches));
        SegmentMUMatches = regexprep(SegmentMUMatches,'[^\-\.0-9]','');
        SegmentMU{j} = str2double(SegmentMUMatches);
        diffSegMU{j} = abs(diff(SegmentMU{j}));
        diffSegMU{j}(end+1) = 0;

        GantryMatches = regexpi(iCOM_data{1},GantryExp,'match','once');
        GantryMatches = GantryMatches(~cellfun(@isempty, GantryMatches));
        GantryMatches = regexprep(GantryMatches,'[^\-\.0-9]','');
        Gantry{j} = str2double(GantryMatches);

        CollimatorMatches = regexpi(iCOM_data{1},CollimatorExp,'match','once');
        CollimatorMatches = CollimatorMatches(~cellfun(@isempty, CollimatorMatches));
        CollimatorMatches = regexprep(CollimatorMatches,'[^\-\.0-9]','');
        Collimator{j} = str2double(CollimatorMatches);

        BeamDescriptionMatches = regexpi(iCOM_data{1},BeamDescriptionExp,'match','once');
        BeamDescriptionMatches = BeamDescriptionMatches(~cellfun(@isempty, BeamDescriptionMatches));
        BeamDescriptionMatches = regexprep(BeamDescriptionMatches,'[^\-\.a-z-A-Z0-9]','');
        BeamDescription{j} = BeamDescriptionMatches;
        
        dMLC_AMatches = regexpi(iCOM_data{1},dMLC_AExp,'match','once');
        dMLC_AMatches = dMLC_AMatches(~cellfun(@isempty, dMLC_AMatches));
        dMLC_AMatches = regexprep(dMLC_AMatches,'[^\-\.0-9]','');
        dMLC_A{j} = str2double(dMLC_AMatches);

        dMLC_BMatches = regexpi(iCOM_data{1},dMLC_BExp,'match','once');
        dMLC_BMatches = dMLC_BMatches(~cellfun(@isempty, dMLC_BMatches));
        dMLC_BMatches = regexprep(dMLC_BMatches,'[^\-\.0-9]','');
        dMLC_B{j} = str2double(dMLC_BMatches);

        dMLC_Matches = regexpi(iCOM_data{1},dMLC_Exp,'match','once');
        dMLC_Matches = dMLC_Matches(~cellfun(@isempty, dMLC_Matches));
        dMLC_Matches = regexprep(dMLC_Matches,'[^\-\.0-9]','');
        dMLC{j} = str2double(dMLC_Matches);
        
        


        k = 1;
        j = 1;
        for i = 1:size(dMLC,2)
            for j = 1:164:length(dMLC{i})
                X1_JAW{i}{k} = dMLC{i}(j,1)*10;
                X2_JAW{i}{k} = dMLC{i}(j+1,1)*10;
                Y1_JAW{i}{k} = dMLC{i}(j+2,1)*10;
                Y2_JAW{i}{k} = dMLC{i}(j+3,1)*10;
                %dMLC_A{i}(j:j+1,1) = '' %get rid of the jaws 
                k = k+1;
                
            end
        k = 1;
        end
        for i = 1:size(dMLC,2)
            for j = 5:164:length(dMLC{i})
                %remove the jaws
                MLC_AB{i}{k} = dMLC{i}(j:j+159,1)*10;
                MLC_A{i}{k} = (MLC_AB{i}{k}(1:2:length(MLC_AB{i}{k}),1));
                MLC_B{i}{k} = (MLC_AB{i}{k}(2:2:length(MLC_AB{i}{k}),1));
                MLC_A{i}{k} = flipud(MLC_A{i}{k});
                MLC_B{i}{k} = flipud(MLC_B{i}{k});
                
                
        
                %dMLC_A{i}(j:j+1,1) = '' %get rid of the jaws 
                k = k+1;
                
            end
        k = 1;
        end
        

      
       
        %now remove the control points that have zero change
        for k = 1: length(diffDelMU)
            ind_cp_mu{k} = find(diffDelMU{k}==0);
            MLC_A{k}(ind_cp_mu{k}) = [ ];
            MLC_B{k}(ind_cp_mu{k}) = [ ];
            X1_JAW{k}(ind_cp_mu{k}) = [ ];
            X2_JAW{k}(ind_cp_mu{k}) = [ ];
            Y1_JAW{k}(ind_cp_mu{k}) = [ ];
            Y2_JAW{k}(ind_cp_mu{k}) = [ ];
            diffSegMU{k}(ind_cp_mu{k}) = [ ];
            diffMU{k}(ind_cp_mu{k}) = [ ];
            diffDelMU{k}(ind_cp_mu{k}) = [ ];
            SegmentMU{k}(ind_cp_mu{k}) = [ ];
            Gantry{k}(ind_cp_mu{k}) = [ ];
            Collimator{k}(ind_cp_mu{k}) = [ ];
            TotalMU{k}(ind_cp_mu{k}) = [ ];
            DeliveryMU{k}(ind_cp_mu{k}) = [ ];

        end
        %now remove the control points that have zero change

       %now we can normalise the mu_del
       for k = 1: length(diffDelMU)
           mu_sum{k} = cumsum(diffDelMU{k});
           mu_max{k} = max(mu_sum{k});
           mu_norm{k} = mu_sum{k}./mu_max{k};
           %mu_sum{k} = [0; diff(mu_sum{k})]
       end
    
    
    %set the parameters for our simulation    
    fluence_grid = ones(400,400);
    fluence_grid_total = zeros(400,400);
    MLC_mask = ones(400,400);
    final_fluence_plan = zeros(400,400);
    %need to set the axes and show the fluce
    

        %setup the MLC parameters
    n_small_leaves = 80;
    leaf_width_large = 5; %mm
    leaf_width_small = 5; %mm
    large_leaf_offset = leaf_width_large;
    leaf_length = 200;
    num_leaves = 80;

for k = 1:length(MLC_A)
    
    %no graphics
    j = 1;
    
    %create the rectangle MLC leaves
    for i = 0:5:395
        A(j) = images.roi.Rectangle('Position',[0,i,200,leaf_width_large],'Color','r','InteractionsAllowed','none');
        A(j).FaceAlpha = 1.0;
        j=j+1;
    end
        
    j = 1;
    %create the rectangle MLC leaves
    for i = 0:5:395
        B(j) = images.roi.Rectangle('Position',[200,i,200,leaf_width_large],'Color','b','InteractionsAllowed','none');
        B(j).FaceAlpha = 1.0;
        j=j+1;
    end
%     %create the rectangle MLC leaves
   
    %create the JAW ROI
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
     %set an original jaw position to return to
    JAW_X1_shift(1) = X1_JAW_ROI(1).Position(1);
    JAW_X2_shift(1) = X2_JAW_ROI(1).Position(1);
    JAW_Y1_shift(1) = Y1_JAW_ROI(1).Position(2);
    JAW_Y2_shift(1) = Y2_JAW_ROI(1).Position(2);
   
    
    
        %for l = 2:2:length(MLC_A{k})-1
        for l = 1:length(MLC_A{k}) 
%               if (diffSegMU{k}(l) < 1) || (diffSegMU{k}(l) > 100)
%                   disp(strcat('found',num2str(k),'/',num2str(l)))
%                   continue
%               end
              
                opening_pos{k}(:,1) = (MLC_A{k}{l}(:));
                opening_pos{k}(:,2) = (MLC_B{k}{l}(:));
                %opening_pos{k} = (fliplr(opening_pos{k}));
               
             
%               figure(100)
%               hold on
%               plot(opening_pos{k}(:,2))

              %these are the max MLC positions used to create the final
              %mask for the most extended leaf
              %these are the max MLC positions used to create the final
              %mask for the most extended leaf
              A_max = min(floor(opening_pos{k}(:,1)));
              B_max = max(ceil(opening_pos{k}(:,2)));
                
              %jaws:=
              opening_pos_x1{k}(1,1) = A_max+2;
              opening_pos_x2{k}(1,1) = B_max-2;
              opening_pos_y1{k}(1,1) = Y2_JAW{k}{l};
              opening_pos_y2{k}(1,1) = -1*Y1_JAW{k}{l};
              
              
              %update the next control point for the jaws
              X1_JAW_ROI(1).Position(1) = X1_JAW_ROI(1).Position(1)+(opening_pos_x1{k}(1,1));
              X2_JAW_ROI(1).Position(1) = X2_JAW_ROI(1).Position(1)+(opening_pos_x2{k}(1,1));
              Y1_JAW_ROI(1).Position(2) = Y1_JAW_ROI(1).Position(2)-(opening_pos_y2{k}(1,1));
              Y2_JAW_ROI(1).Position(2) = Y2_JAW_ROI(1).Position(2)-(opening_pos_y1{k}(1,1));
           
               for m = 1:length(MLC_A{k}{l})% MLC positions for one control point control point
    
                      
                    %update current position to next position
                    A(m).Position(1) = A(m).Position(1)+(opening_pos{k}(m,1)); %leaf offset A
                    B(m).Position(1) = B(m).Position(1)+(opening_pos{k}(m,2)); %leaf offset A
                   

                     mask_A = ~createMask(A(m),fluence_grid);
                     mask_B = ~createMask(B(m),fluence_grid);
    
                     MLC_mask = MLC_mask.*(mask_A.*mask_B);
                     
                
               %multiWaitbar( 'Processing MLC Pair...', (m/length(MLC_A{k}{l})), 'Color', 'g');
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
            
            fluence_grid_total = (fluence_grid_total+(fluence_grid.*JAW_MLC_mask{l}.*diffDelMU{k}(l)));
            %fluence_grid_total(fluence_grid_total~=0) = (fluence_grid_total(fluence_grid_total~=0)+diffDelMU{k}(l));
            %fluence_grid_total = fluence_grid_total - 1;
            %fluence_grid_total(fluence_grid_total<=0) = 0;
%             if abs(Gantry{k}(l)) < 90
%                 disp('flippin g')
%                 fluence_grid_total = fliplr(fluence_grid_total);
%             end

            %fluence_grid_total(fluence_grid_total~=0) = (fluence_grid_total(fluence_grid_total~=0)-1);
            %fluence_grid_total(fluence_grid_total~=0) = (fluence_grid_total(fluence_grid_total~=0)-1);
            %fluence_grid_total(fluence_grid_total~=0) = (fluence_grid_total(fluence_grid_total~=0);
            for i = 1:length(A)
                A(i).Position(1) = MLC_A_shift(i);
                B(i).Position(1) = MLC_B_shift(i);
            end
            
            %reset the JAW ROI
            X1_JAW_ROI(1).Position(1) = JAW_X1_shift(1);
            X2_JAW_ROI(1).Position(1) = JAW_X2_shift(1);
            Y1_JAW_ROI(1).Position(2) = JAW_Y1_shift(1);
            Y2_JAW_ROI(1).Position(2) = JAW_Y2_shift(1);
            fluence_grid = ones(400,400); 
            %drawnow
       
        
       
        if Collimator{k}(l)<180.1 %clockwise rotation
            %final_fluence_beam{k} = flipud(imrotate(fluence_grid_total,1*CollAngle(k,1),'bilinear','crop'));
            final_fluence_beam{k} = (imrotate(fluence_grid_total,1*Collimator{k}(l),'bilinear','crop'));
            %final_fluence_plan = imrotate(final_fluence_plan,-1*CollAngle(k,1),'bilinear','crop');
        else %counter clockwise rotation
            %final_fluence_beam{k} = flipud(imrotate(fluence_grid_total,-1*CollAngle(k,1),'bilinear','crop'));
            final_fluence_beam{k} = (imrotate(fluence_grid_total,-1*Collimator{k}(l),'bilinear','crop'));
            %final_fluence_plan = imrotate(final_fluence_plan,1*CollAngle(k,1),'bilinear','crop');
        end
        
        

    %add them to the final plan fluence and then reset the coll
        
        end
    fluence_grid_total = zeros(400,400);   
    %final_fluence_beam{k} = (fluence_grid_total);
    
    final_fluence_plan = final_fluence_plan + final_fluence_beam{k};
    %final_fluence_plan = medfilt2(final_fluence_plan,[4 4]);
    
    IDf = ID{1}{1}(1,4:end)
    Rx{1}{1}
    Rxf = Rx{1}{1}
    Beamf = Beam{k}{1}
    Machinef = Machine{1}{1}(1,4:end)
    fileNameLog = strcat(OutFolder,IDf,'_',Rxf,'_',Machinef,'_',fName,'_TOTAL','_DEL','.csv') %this is the filename for the Watcher / Matcher
    %writematrix(final_fluence_beam{k},strcat(OutFolder,IDf,'_',Rxf,'_Beam_',Beamf,'_',Machinef,'_',fName,'_DEL','.csv'));
  
    
end

    writematrix(final_fluence_plan,strcat(OutFolder,IDf,'_',Rxf,'_',Beamf,'_',Machinef,'_',fName,'_TOTAL','_DEL','.csv'))


        

     
      

    %app.UITable.Data = [TableData; {datestr(datetime('now')),IDf, Rxf, Beamf, Machinef, fileNameLog}]
    
    %call the matcher and kick off the matching process
    %MatcherLogToPlan(app,IDf,Rxf,'TOTAL',Machinef,fileNameLog,OutFolder,1,4,10)
end   
