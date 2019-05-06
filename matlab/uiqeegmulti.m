function [sf] = uiqeegmulti()

%%
sf=[];

%% Get input files
[file, path] = uigetfile('*.easy', 'Select Pre,Inter,Pose easy files', 'MultiSelect', 'on');
if length(file) ~= 3
    error('Expecting three input files.');
end

% Labels used in legends for spectra, band plots. These should be in the
% same order as the order in which the files are selected. 
myLegend = {'Pre', 'Inter', 'Post'};

% Channel numbers to analyze. Should be in the same order as the
% order in which the files are chosen.
myChannels = [7 7 7];

sf = qeegmulti({ [path file{1}], [path file{2}], [path, file{3}] }, ...
                'Channels', myChannels, 'Legend', myLegend);

return
