function MatcherLogToPlan(app,LogID,LogRx,LogBeam,LogMachine, inputFile, OutFolder, dta, dd, thresh)
%function [outputArg1,outputArg2] = Matcher(delcsvpath)
%MATCHER Summary of this function goes here
%   Detailed explanation goes here

%1. look in the ouput directory for files that match the plan filename

%get the list of files in the dir
filesList = dir(OutFolder)
sstr1=strsplit(inputFile,'_');
patientID = extractBefore(inputFile,'_');
patientID = regexprep(patientID,'[a-zA-Z\s \: \\]','');
disp(patientID)
%tf = strcmp(filesList.name(1),inputFile)
j = 1
for i = 1:length(filesList)

    
    TF_plan = contains(filesList(i).name,'_TOTAL_PLAN');
    
    if TF_plan == 1
        disp('Found a plan file')
        %get the matching fileName
        TF_planID = contains(filesList(i).name,patientID,'IgnoreCase',true);
                if TF_planID == 1
                    fileNameforComp{j} = filesList(i).name
                    filePathforComp{j} = filesList(i).folder
                    filePathRef = strcat(filePathforComp{j},'\',fileNameforComp{j})
                    disp(fileNameforComp{j});

                    %parfeval(app.p,@GammaEval,0,app, app.PlanID, app.PlanRx, app.PlanBeam, app.PlanMachine, inputFile, filePathRef, dta, dd, thresh)
                    GammaEval(app,LogID,LogRx,LogBeam,LogMachine, inputFile, filePathRef, dta, dd, thresh)
                    j = j+1
                end
       
            
        
    end

    
end
end

