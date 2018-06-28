%% this script filters and thresholds stimulation data
    clear
    pwd = cd;

    inputData.folderpath= 'C:\Users\jts3256\Desktop\20161104\'; % must have \ at the end
%     inputData.folderpath = 'C:\Users\ttp603\Desktop\data\sorting\butter\Butter_20180504_cuneate\';
%     inputData.folderpath = 'D:\Lab\Data\StimArtifact\testData\';
%     inputData.mapFile='mapFileY:\limblab\lab_folder\Animal-Miscellany\Butter_17D2\Right Cuneate Implant 2018_03_06\SN 6250-001799.cmp';
%     inputData.mapFile = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
    inputData.mapFile = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Chips_12H1\map_files\left S1\SN 6251-001455.cmp';
    inputData.task='taskRT';
    inputData.ranBy='ranByJoseph'; 
    inputData.array1='arrayLeftS1'; 
    inputData.monkey='monkeyChips';
    inputData.labnum = 6;
    
    inputData.dukeBoardChannel = -1;
    inputData.dukeBoardLabel = 'ainp15';

    inputData.badChList=0;
    inputData.interpulse=.000053;%in s
    inputData.pWidth1=.0002;
    inputData.pWidth2=.0002;

    % functionName='processStimArtifactData';

    inputData.artifactDataTime = 10; % in ms
    inputData.presample = 100;

    inputData.preOffset = 22;
    inputData.postOffset = 25;

    inputData.moreThanOnePulsePerWave = 0;
    inputData.numPulses = 10;
    inputData.pulseFrequency = 100;

    inputData.thresholdMult = 3.5;
    inputData.artifactSkip = 1;
    inputData.maxAmplitude = 1000; % in uV

    inputData.maxChunkLength = 5000*30; % 5 second chunk maximum
%% get spike crossing for all files, write a nev file
    cd(inputData.folderpath)
    fileList = dirSorted('*.ns5');
    outputData = [];
    % save input data
    save(strcat(fileList(1).name(1:end-4),'_inputData.mat'),'inputData');

    nevDataAll = [];
    durationAll = 0;
    % process data
    for f = 1:numel(fileList)
        warning('off')
        inputData.filename = fileList(f).name;

        [~,outputData] = filterAndThresholdData(inputData);
        stimInfo = outputData.stimInfo;
        stimInfo.chanSent = outputData.waveforms.chanSent';
        stimInfo.waveSent = outputData.waveforms.waveSent';
        stimInfo.parameters = outputData.waveforms.parameters;
       
        save([inputData.folderpath,inputData.filename(1:end-4),'_outputData.mat'],'outputData','-v7.3');
        save([inputData.folderpath,inputData.filename(1:end-4),'_stimInfo.mat'],'stimInfo','-v7.3');
        % append nev data
        if(f == 1)
            nevDataAll = outputData.nevData;
        else
            numAdd = numel(outputData.nevData.ts);
            nevDataAll.ts(end+1:end+numAdd,:) = outputData.nevData.ts + durationAll;
            nevDataAll.waveforms(end+1:end+numAdd,:) = outputData.nevData.waveforms;
            nevDataAll.elec(end+1:end+numAdd,:) = outputData.nevData.elec;
        end
        durationAll = durationAll + outputData.duration;
    end

    
    
    % save('','artifactData')

    % write nev file
    disp('writing nev file')

    packetWidth = 104;
    underscoreIdx = strfind(fileList(1).name,'_');
    filename = strcat(fileList(1).name(1:underscoreIdx(2)),'_merged');
    mapFilename = inputData.mapFile(8:end);
    comments = '';
    writeNEV(nevDataAll, packetWidth, filename, mapFilename, comments )

    cd(pwd);
    disp('DONE -- CAN CONTINUE')
    warning('on')

%% sort *_merged and call it *_merged-s, rename ns1 as well (?)

%% load in *_merged-s, split apart the units,
%% save the units to a new nev file as well as to output data and a stimInfo file
% the new nev file (*_spikesExtracted*), the new ns5 (*_spikesExtracted*), and the stimInfo mat file are needed for
% further processing. Everything else is excess

    disp('started saving unit files')

    pwd = cd;
    cd(inputData.folderpath);
    outputDataFileList = dirSorted('*outputData.mat');
    nevFileList = dirSorted('*_merged-s.NEV*');
    NEVname = nevFileList(1).name; % grab the first one
    
    NEV_dataAll = openNEV('read', [inputData.folderpath NEVname],'nosave');
    durationAll = 0;
    totalUnits = 0;
    for f = 1:numel(outputDataFileList)
        load(outputDataFileList(f).name);
        units = [];
        unitsMask = [];
        
        % split back into individual files

        unitsMask = double(NEV_dataAll.Data.Spikes.TimeStamp)/30000 - durationAll < outputData.duration & ...
            double(NEV_dataAll.Data.Spikes.TimeStamp)/30000 - durationAll > 0;
        units.ts = double(NEV_dataAll.Data.Spikes.TimeStamp(unitsMask))/30000 - durationAll;
        units.elec = NEV_dataAll.Data.Spikes.Electrode(unitsMask);
        units.label = NEV_dataAll.Data.Spikes.Unit(unitsMask);
        units.waveform = NEV_dataAll.Data.Spikes.Waveform(:,unitsMask)*0.254;
        
        % undo any duration adding do to resets
        for resetIdx = 1:numel(outputData.DataDurationSec)-1 % to prevent an extra reset
            mask = units.ts > outputData.DataDurationSec(resetIdx);
            units.ts(mask) = units.ts(mask) - outputData.DataDurationSec(resetIdx);
            % stim on info is already adjusted (incorrectly, ignore this)
        end
        % load normal nev
        NEV_dataSingle = openNEV('read',[inputData.folderpath outputDataFileList(f).name(1:end-15) '.nev'],'nosave');
        
        % replace spike data in the normal NEV
        NEV_dataSingle.Data.Spikes.TimeStamp = uint32(units.ts*30000); %uint32
        NEV_dataSingle.Data.Spikes.Electrode = units.elec;
        NEV_dataSingle.Data.Spikes.Waveform = units.waveform;
        NEV_dataSingle.Data.Spikes.Unit = units.label;
%         NEV_dataSingle.DataDuration = outputData.DataPoints;
%         NEV_dataSingle.DataDurationSec = outputData.DataDurationSec;
        
        % save normal nev with a new name
        saveNEV(NEV_dataSingle,[inputData.folderpath outputDataFileList(f).name(1:end-15), '_spikesExtracted.nev'],'noreport');
        
        
        % update duration for splitting files
        durationAll = durationAll + outputData.duration;
    end
    
    disp('done replacing spike info')