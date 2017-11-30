%% Volume occupancy of one mask within another
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
% *Loads two image stacks containing binari masks: the content mask contains
% a signal widely distributed in the volume. The container mask defines a
% restricted zone in which quantify the content. Then calculate the 
% intersection mask of those two and report occupancy statistics.*
%
% For instance, the content mask could be a synaptic labelling widespread 
% within a neuropil layer and the container mask could define the regions 
% in which axons of the neurons of interest stratity.
% This program will calculate how many of the synpses belong to that 
% particular neuron type labeled in the container mask.
%
% The program does the following main operations:
%
% # Ask user to load the container mask
% # Ask user to load the content mask
% # Ask user to set resolution and whether to save intersetion mask or not
% # Calculate the content-within-container mask (and save it if desidered)
% # Display on-screen absolute volumes and percentage occupancy
%
% *Input:* 
%
% * Container mask (i.e. the mask of neuron's axon terminal)
% * Content mask (i.e. the mask of puncta in the entire synaptic layer)
%
% *Output :* 
%
% * Content within container mask (.tif file)
% * Volume occupancy statistics
%          (i.e. mask of puncta within the neuron's axon terminal)
%
% *Dependencies:*
%
% * textprogressbar.m (a fast implementation of a console progress bar)


% Load the container mask
[FileName, PathName] = uigetfile('*.tif', 'Load the mask of the container signal');
tmpContainerImInfo = imfinfo([PathName FileName]);
tmpContainerImSize = [tmpContainerImInfo(1).Height tmpContainerImInfo(1).Width length(tmpContainerImInfo)];
textprogressbar(['Loading container image stack: "' FileName '" ...']);
tmpContainerStack = zeros(tmpContainerImSize(1), tmpContainerImSize(2), tmpContainerImSize(3));
for i = 1:tmpContainerImSize(3)
    tmpContainerStack(:,:,i)=imread([PathName FileName], i);
    textprogressbar(100*i/tmpContainerImSize(3)); % update progress bar
end
tmpContainerStack = tmpContainerStack / max(max(max(tmpContainerStack))); % Normalize mask values to zeros and ones
textprogressbar('DONE');

% Load the content mask
[FileName, PathName] = uigetfile('*.tif', 'Load the mask of the content signal');
tmpContentImInfo = imfinfo([PathName FileName]);
tmpContentImSize = [tmpContentImInfo(1).Height tmpContentImInfo(1).Width length(tmpContentImInfo)];
textprogressbar(['Loading content image stack: "' FileName '" ...']);
tmpContentStack = zeros(tmpContentImSize(1), tmpContentImSize(2), tmpContentImSize(3));
for i = 1:tmpContentImSize(3)
    tmpContentStack(:,:,i)=imread([PathName FileName], i);
    textprogressbar(100*i/tmpContentImSize(3)); % update progress bar
end
tmpContentStack = tmpContentStack / max(max(max(tmpContentStack))); % Normalize mask values to zeros and ones
textprogressbar('DONE');

% Input user-defined parameters
tmpPrompt = {'Specify X-Y resolution of the image stacks (µm/pixel):',...
             'Specify Z resolution of the image stacks (µm/pixel):',...
             'Save mask of content within container? 0=no, 1=yes'};
tmpAns = inputdlg(tmpPrompt, 'Image resolution',[1 40],{'0.097','0.3','1'});

tmpResXY    = str2double(tmpAns{1});   % number of gamma-fit standard devs.
tmpResZ     = str2double(tmpAns{2});   % fit mode (z-depth independency)
tmpSaveMask = str2double(tmpAns{3});   % save the intersection mask? 

% Calculate statistics of compared masks

% intersect the two maks to create the mask of content within container
tmpContentWithinContainerStack = tmpContentStack.*tmpContainerStack;

% save intersection mask as tif if we are in debug mode
if tmpSaveMask == 1
    [FileName,PathName,~] = uiputfile('*.tif','Name&Save the tif file.');
    textprogressbar(['Saving Content within Container Mask in "' FileName '": ']);
    imwrite(tmpContentWithinContainerStack(:,:,1), [PathName FileName], 'tif', 'compression', 'lzw'); %write first z plane
    if size(tmpContentWithinContainerStack,3) > 1 %write the rest of the z planes
        for i=2:size(tmpContentWithinContainerStack,3)
            imwrite(tmpContentWithinContainerStack(:,:,i), [PathName FileName], 'tif', 'compression', 'lzw', 'WriteMode', 'append');
            textprogressbar(100*i/size(tmpContentWithinContainerStack,3)); % update progress bar
        end
    end
    textprogressbar('DONE');
end

tmpVoxelVolume = tmpResXY*tmpResXY*tmpResZ;

% Print statistics on screen
fprintf('\nStatistics of the processed images:');
fprintf(['\nVolume of compared image stacks    (µm^3):' num2str(numel(tmpContainerStack)*tmpVoxelVolume)]);
fprintf(['\nVolume of the container mask       (µm^3):' num2str(numel(find(tmpContainerStack))*tmpVoxelVolume)]);
fprintf(['\nVolume of the entire content mask  (µm^3):' num2str(numel(find(tmpContentStack))*tmpVoxelVolume)]);
fprintf(['\nVolume of content within container (µm^3):' num2str(numel(find(tmpContentWithinContainerStack))*tmpVoxelVolume)]);
fprintf(['\nPercent of container occupied by content :' num2str(100*numel(find(tmpContentWithinContainerStack))/numel(find(tmpContainerStack)))]);
fprintf('\n');

% Cleanings
clear tmp* ans FileName PathName FilterIndex i;

%% Change log
% _*Version 2.4*             created on 2017-11-30 by Luca Della Santina_ 
%
%  + Added copyright and license information
%  % Resolved all MATLAB warnings
%
% _*Version 2.2*             created on 2017-09-08 by Luca Della Santina_ 
%
%  % Fixed wrong calculation of absolute volumes (using Xres*Yres*Zres)
%
% _*Version 2.1.1*             created on 2017-09-08 by Luca Della Santina_ 
%
%  + Reformatted documentation using proper MATLAB markup format 
%
% _*Version 2.1*               created on 2017-04-15 by Luca Della Santina_ 
%
%  + Report the absolute volume of the content within the container mask
%
% _*Version 2.0*               created on 2017-04-14 by Luca Della Santina_ 
%
%  + Computes absolute volumes and percent volume occupancy
%  + Save intersection mask of content within container
%  + Allows user to use custom XYZ voxel resolution
%
% _*Version 1.0*               created on 2017-04-01 by Luca Della Santina_ 