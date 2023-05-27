function UpdateTable(app,LogID,LogRx,LogBeam,LogMachine,numpass,filePathRef,filePathEval)
% Compute statistics
LogID
LogBeam
LogMachine
numpass

TableData = app.UITable.Data
%current data
%ADate, ID, plan, Beam, GPR, Plan path, Log Path
app.UITable.Data = [TableData; {datestr(datetime('now')), LogID, LogRx, LogBeam, LogMachine, numpass*100, filePathRef, filePathEval}]


%UPDATETABLE Summary of this function goes here
%   Detailed explanation goes here

end

