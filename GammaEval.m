function [numpass, avg, GammaMap, image1, image2, filePathEval, filePathRef, dta, dd, thresh] = GammaEval(app, LogID, LogRx, LogBeam, LogMachine, filePathEval, filePathRef,  dta, dd, thresh)
%GAMMEVAL Summary of this function goes hereGammaEval(app, ID, Plan, Beam, Machine, inputFile, filePathRef, dta, dd, thresh)
%   Detailed explanation goes here
filePathRef
filePathEval
LogID
LogBeam
LogMachine
image2 = readmatrix(filePathEval,'FileType','text','Delimiter',','); %delivery
image1 = readmatrix(filePathRef,'FileType','text','Delimiter',','); %plan

%normalise
image1 = image1./(max(max(image1)));
image2 = image2./(max(max(image2)));

debuglevel = 2;
res_x = 1;
res_y = 1;
DTA_tol = dta;
radlim = DTA_tol * 1.5;
rad = min(ceil(radlim/res_x),ceil(radlim/res_y));
MaxVal = max(image1(:));
FE_thresh = thresh/100;
Dose_tol = dd/100;
%MaxVal = max(max(Image1));
%[Im1_ydim, Im1_xdim] = size(Image1);
Mask = zeros(size(image1));
crit_val = FE_thresh*MaxVal;

%[Im1_ydim Im1_xdim] = size(Image1);

% Use the Resampled image for this - EPID no spikes?
Mask(image2>crit_val) = 1;

Dose_tol = Dose_tol*MaxVal; % GLOBAL - PERCENTAGE OF MAX DOSE


if debuglevel > 1
    fprintf('Maximum dose for Gamma: %2.1f MU\n',MaxVal);
    fprintf('Dose Tolerance for Gamma: %2.1f MU\n',Dose_tol);
    fprintf('DTA: %2.1f\n', DTA_tol);
    fprintf('FE_thresh: %2.1f\n', FE_thresh);
    fprintf('rad: %2.1f\n', rad);
    fprintf('Image 1 val: %2.1f\n',image1(round(size(image1,1)/2), round(size(image1,2)/2)));
    fprintf('\n');
end

% VECTORIZED CALCULATION STARTS HERE

% GammaMapsub will carry the calculated gamma values for the truncated
% images. GammaMap2 will be the Gamma values for the full image.
GammaMapsub = NaN;
GammaMap = zeros(size(image1)).*image1;  % casts to gpu if present

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
if max_x > size(image1,2)
    max_x = size(image1,2);
end
if max_y > size(image1,1)
    max_y = size(image1,1);
end
% Truncate the images to avoid needless calculations
Im1 = image1(min_y:max_y,min_x:max_x);
Im2 = image2(min_y:max_y,min_x:max_x);
% Shift the image by varying amounts. Determine the minimum gamma value
% for all shifts
for i=-rad:rad
    for j=-rad:rad
        % circshift function wraps elements from top to bottom as necessary
        % The entire image is shifted at once
        Im2_shift = circshift(Im2,[i j]);
        dist = sqrt((res_y*i)^2 + (res_x*j)^2);
        DoseDiff = Im2_shift - Im1;
        % Compute the gamma map for this particular shift value
        Gamma_temp = sqrt((dist./DTA_tol).^2 + (DoseDiff./Dose_tol).^2);
        % Accumulate the map of the minimum values of gamma at each point
        GammaMapsub = min(GammaMapsub,Gamma_temp);
    end
end
% Put the truncated gamma map back into its proper location within the full
% gamma map
GammaMap(min_y:max_y,min_x:max_x) = GammaMapsub;
% Remove any edge effects from the circular shifting by multiplying by the
% mask values. This will negate any calculated gamma values around the
% edges of the distribution where this effect would arise
GammaMapNoMask = GammaMap;
GammaMap = GammaMap.* Mask;
% Ensure that NaN values outside the mask do not affect the calculation
GammaMap(~Mask) = 0.0;

% Compute statistics
numWithinField = nnz(Mask)
numpass = nnz(GammaMap<1 & Mask)./numWithinField
app.numpass = numpass;
avg = sum(GammaMap(:))./numWithinField
app.avg = avg;

TableData = app.UITable.Data
%current data
%ADate, ID, plan, Beam, GPR, Plan path, Log Path
app.UITable.Data = [TableData; {datestr(datetime('now')), LogID, LogRx, LogBeam, LogMachine, app.numpass*100, filePathRef, filePathEval}]



end

