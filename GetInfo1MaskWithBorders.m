%% Binary Mask statistics (single mask)
% |Copyright 2017, Luca Della Santina|
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% This software is released under the terms of the GPL v3 software license
%
% *Loads a binary mask and display its properties
%
% *Input:* 
%
% * Mask file (i.e. the mask of neuron's axon terminal)
% * 2 sets of XYZ coordinates (i.e. defining IPL boundaries)
%
% *Output :* 
%
% * Absolute volume of image
% * Absolute volume of mask
% * Percent occupancy of the image by mask
% * Percent stratification of the mask within borders (min-max)
%
% *Dependencies:*
%
% * textprogressbar.m (a fast implementation of a console progress bar)
% * gridfit.m (for fitting a surface to a set of XYZ coordinates)

% Load the container mask
disp('---- GetMaskInfo (single mask with borders) ver. 2.4 ----');
[FileName, PathName] = uigetfile('*.tif', 'Load the mask of the container signal');
tmpMaskImInfo = imfinfo([PathName FileName]);
tmpMaskImSize = [tmpMaskImInfo(1).Height tmpMaskImInfo(1).Width length(tmpMaskImInfo)];
textprogressbar(['Loading image stack: "' FileName '" ...']);
tmpMaskStack = zeros(tmpMaskImSize(1), tmpMaskImSize(2), tmpMaskImSize(3));
for i = 1:tmpMaskImSize(3)
    tmpMaskStack(:,:,i)=imread([PathName FileName], i);
    textprogressbar(100*i/tmpMaskImSize(3)); % update progress bar
end
tmpMaskStack = tmpMaskStack / max(max(max(tmpMaskStack))); % Normalize mask values to zeros and ones
textprogressbar('DONE');

% Get dimensions of the first image (assumes each channel stored as individual 3D .tif file coming from same acquisition)
tmpImSizeZ = numel(tmpMaskImInfo);
tmpImSizeX = tmpMaskImInfo.Width;
tmpImSizeY = tmpMaskImInfo.Height;

% Retrieve XY and Z resolution from TIFF image descriptor or ask user if not available
if isempty(tmpMaskImInfo(1).XResolution)
    % ask user to input image resolution because not store in the image
    tmpPrompt = {'X-Y resolution (µm/pixel):', 'Z resolution (µm/pixel):','Offset from INL points (µm)', 'Offset from GCL points (µm)'};
    tmpAns = inputdlg(tmpPrompt, 'Image resolution',[1 40],{'0.189','0.3','8','10'});
    tmpResXY = str2double(tmpAns{1}); 
    tmpResZ = str2double(tmpAns{2});
    offset1 = str2double(tmpAns{3});
    offset2 = str2double(tmpAns{4});
    offset1 = ceil(offset1/tmpResZ); % convert offset microns to pixels
    offset2 = ceil(offset2/tmpResZ); % convert offset microns to pixels
else
    tmpResXY = num2str(1/tmpMaskImInfo(1).XResolution);
    if contains(tmpMaskImInfo(1).ImageDescription, 'spacing=')
        tmpPos = strfind(tmpMaskImInfo(1).ImageDescription,'spacing=');
        tmpResZ = tmpMaskImInfo(1).ImageDescription(tmpPos+8:end);
        tmpResZ = regexp(tmpResZ,'\n','split');
        tmpResZ = tmpResZ{1};
    else
        tmpResZ = '0.3'; % otherwise use default value
    end
end
tmpVoxelVolume = tmpResXY*tmpResXY*tmpResZ; % store the real size of voxels


% Load the points delimiting the INL/IPL surface, then surfaceFit them

fitDownScaling = 2; % Downscale 2x the X and Y coordinates for efficient surfacefit
if exist('Spots1.mat', 'file')
   load('Spots1.mat'); 
else
    [FileNamePoints, PathNamePoints] = uigetfile('*.mat', 'Select INL limit points', 'MultiSelect', 'off');
    load([PathNamePoints FileNamePoints]);
end
SpotsXYZ=double(SpotsXYZ);
x=SpotsXYZ(:,1);
y=SpotsXYZ(:,2);
z=SpotsXYZ(:,3);
gx=1:fitDownScaling:tmpImSizeX;
gy=1:fitDownScaling:tmpImSizeY;
fprintf('Fitting INL limit points with a surface...');
tic;
limitsINL=gridfit(x,y,z,gx,gy, 'smoothness', 75); % default smoothness=1
limitsINL=ceil(limitsINL)+offset1;
fprintf(['DONE in ' num2str(toc) ' seconds \n']);

% Load the points delimiting the IPL/GCL surface, then surfaceFit them

if exist('Spots2.mat', 'file')
    load('Spots2.mat'); 
else
    [FileNamePoints, PathNamePoints] = uigetfile('*.mat', 'Select GCL limit points', 'MultiSelect', 'off');
    load([PathNamePoints FileNamePoints]);
end
SpotsXYZ=double(SpotsXYZ);
x=SpotsXYZ(:,1);
y=SpotsXYZ(:,2);
z=SpotsXYZ(:,3);
gx=1:fitDownScaling:tmpImSizeX;
gy=1:fitDownScaling:tmpImSizeY;
fprintf('Fitting GCL limit points with a surface...');
tic;
limitsGCL=gridfit(x,y,z,gx,gy, 'smoothness', 75); % default smoothness=1
limitsGCL=ceil(limitsGCL)-offset2;
fprintf(['DONE in ' num2str(toc) ' seconds \n']);

% Convert the mask into a matrix of IPL stratification depth (per voxel)

tmpMaskIdx = find(tmpMaskStack)';
tmpMaskIdxSize= numel(tmpMaskIdx);
fprintf('Calculating IPL stratification level... ');
tic;
for i=tmpMaskIdx
    [x, y, z] = ind2sub([tmpImSizeX, tmpImSizeY, tmpImSizeZ], i);
    tmpLimitINL = limitsINL(ceil(x/fitDownScaling), ceil(y/fitDownScaling));
    tmpLimitGCL = limitsGCL(ceil(x/fitDownScaling), ceil(y/fitDownScaling));
    tmpZperc = 100*(z-tmpLimitINL)/(tmpLimitGCL-tmpLimitINL);
    tmpMaskStack(i) = tmpZperc;
end
fprintf(['DONE in ' num2str(toc) ' seconds \n']);

% Bin stratifications by IPL percent (stepping = 1 percent)
tmpDensity=zeros(1, 100);
fprintf('Count masked voxel along IPL depth... ');
tic;
for i = 1:100
    tmpDensity(i) = sum(sum(sum(tmpMaskStack >(i-1) & tmpMaskStack <=i)));
    %Old version (not optimized because find is a slower function than sum) 
    %tmpDensity(i) = numel(find(tmpMaskStack >(i-1) & tmpMaskStack <=i));
end
tmpDensityNorm=tmpDensity/max(tmpDensity);
fprintf(['DONE in ' num2str(toc) ' seconds\n']);

% Plot dot density distribution as a function of Volume depth.
tmpH = figure('Name', 'Binned distribution along Z');
set(tmpH, 'Position', [100 200 1200 500]);
set(gcf, 'DefaultAxesFontName', 'Arial', 'DefaultAxesFontSize', 12);
set(gcf, 'DefaultTextFontName', 'Arial', 'DefaultTextFontSize', 12);

subplot(1,2,1);
hold on;
tmpY = tmpDensity;
tmpX = 1:100;
plot(tmpX, tmpY, 'k', 'MarkerSize', 8);

box off;
set(gca, 'color', 'none',  'TickDir','out');
ylabel('Number of voxels');
xlabel('Volume depth percentage');

subplot(1,2,2);
hold on;
tmpY = tmpDensityNorm;
tmpX = 1:100;
plot(tmpX, tmpY, 'k', 'MarkerSize', 8);

box off;
set(gca, 'color', 'none',  'TickDir','out');
ylabel('Relative voxel density');
xlabel('Volume depth percentage');

[~,FileName,~] = fileparts(FileName);
xlswrite([FileName 'StratNormalized.xls'],[1:100; tmpDensityNorm]');
xlswrite([FileName 'StratAbsolute.xls'],[1:100; tmpDensity]');

% Print statistics on screen
disp(' ');
disp(['statistics of mask "' FileName '"']);
disp(['Image absolute volume       (µm^3): ' num2str(ceil(numel(tmpMaskStack)*tmpVoxelVolume))]);
disp(['Mask absolute volume        (µm^3): ' num2str(ceil(numel(find(tmpMaskStack))*tmpVoxelVolume))]);
disp(['Mask volume occupancy (% of image): ' num2str(100* numel(find(tmpMaskStack))/numel(tmpMaskStack))]);
disp(['Mask stratification level  (% IPL): ' num2str(find(tmpDensity, 1 )) ' - ' num2str(find(tmpDensity, 1, 'last'))]);

% Cleanings
clear tmp* ans FileName PathName fitDownScaling gx gy i limits* x y z SpotsXYZ ext offset*;

%% Change log
% _*Version 2.4*               created on 2017-11-30 by Luca Della Santina_ 
%
%  + Speed optimization in calculating percentages (20x)
%  + Speed optimization in calculating densities (1.5x)
%  % Removed all MATLAB warnings
%
% _*Version 2.3.2*             created on 2017-11-28 by Luca Della Santina_ 
%
%  + User can add custom spacing from the points to delimit IPL region
%  % Changed default image resolution to 0.189x0.189x0.3 µm
%  % Fixed bug in calculation of IPL stratification min/max
%  % Add image FileName as prefix to output xls result files