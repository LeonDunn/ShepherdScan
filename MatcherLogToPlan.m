function [filePathRefm, filePathEvalm] = MatcherLogToPlan(app,LogIDf,LogRxf,LogBeamf,LogMachinef, inputFile, OutFolder)
%function [outputArg1,outputArg2] = Matcher(delcsvpath)
%MATCHER Summary of this function goes here
%   Detailed explanation goes here

%1. look in the ouput directory for files that match the plan filename
LogIDf
LogRxf
LogBeamf
LogMachinef

%get the list of files in the dir
filesList = dir(OutFolder);
sstr1=strsplit(inputFile,'_');
patientID = extractBefore(inputFile,'_');
patientID = regexprep(patientID,'[a-zA-Z\s \: \\]','');

%tf = strcmp(filesList.name(1),inputFile)
j = 1;
    for i = 1:length(filesList)

    
        TF_plan = contains(filesList(i).name,'_TOTAL_PLAN');
        
        if TF_plan == 1
            disp('Found a plan file')
            %get the matching fileName
            TF_planID = contains(filesList(i).name,patientID,'IgnoreCase',true);
                    if TF_planID == 1
                            fileNameforComp{j} = filesList(i).name;
                            filePathforComp{j} = filesList(i).folder;
                            filePathRefm = strcat(filePathforComp{j},'\',fileNameforComp{j});
                            disp(fileNameforComp{j})
                            %fcn_gamma = parfeval(app.p,@GammaEval,1,app, LogID, LogRx, LogBeam, LogMachine, inputFile, filePathRef, dta, dd, thresh);
                            filePathEvalm = inputFile;
                            app.filePathRef = filePathRefm;
                            app.filePathEval = filePathEvalm;
                            j = j+1;
                   end
    
        end

    
    end
end

