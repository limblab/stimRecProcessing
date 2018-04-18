%% process stimulation artifacts:
% this script will grab the artifact related data on all channels and will
% grab an acausal filtered version of the artifact as well. 

% this script does not do any thresholding of the data.

pwd = cd;
folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\Han_20180310_chic201802\';
functionName='processStimArtifact_filter'; % includes acausal filtered data

%%
inputData.saveFigures = 0; % do we want to save the figures we make
inputData.makePlots = 0; % do we want to make figures

% if using the duke board, specify the channel it was connected to and the
% name of the analog channel the data is on
inputData.dukeBoardChannel = -1;
inputData.dukeBoardLabel = 'ainp15';

% input data to make a cds
inputData.task='taskCObump';
inputData.ranBy='ranByTucker'; 
inputData.array1='arrayS1'; 
inputData.monkey='monkeyHan';
inputData.mapFile='mapFileR:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp'; % mapfile location

% channel and stim data
inputData.badChList=[];
inputData.interpulse=.000053;%in s
inputData.pWidth1=.0002;
inputData.pWidth2=.0002;

% plotting data
inputData.windowSize=30*10;%in points
inputData.presample=100;%in points
inputData.plotRange=8;%in mV
inputData.plotRangeFiltered=0.5;%in mV
inputData.lab=6;
inputData.useSyncLabel=[];

% run the function
[output_data] = runDataProcessing(functionName,folderpath,inputData);
disp('done processing')
cd(pwd);


