%% this script filters and thresholds stimulation data
    clear
    pwd = cd;
    inputData.folderpath= 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\testingCode\'; % must have \ at the end
%     inputData.folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\testingCode\';
%     inputData.folderpath = 'D:\Lab\Data\StimArtifact\testData\';
    inputData.mapFile='mapFileR:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
    % inputData.mapFile = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Chips_12H1\map_files\left S1\SN 6251-001455.cmp';

    inputData.task='taskCObump';
    inputData.ranBy='ranByJoseph'; 
    inputData.array1='arrayLeftS1'; 
    inputData.monkey='monkeyHan';
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
%% generates _cds and _nevData files, also writes nev file
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
        save([inputData.folderpath,inputData.filename(1:end-4),'_outputData.mat'],'outputData');
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

%% sort *_merged and call it *_merged-s

%% load in *_merged-s, save the units to each nev file as well as to output data
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
        for resetIdx = 1:numel(outputData.DataDurationSec)
            mask = units.ts > outputData.DataDurationSec(resetIdx);
            units.ts(mask) = units.ts(mask) - outputData.DataDurationSec(resetIdx);
            % stim on info is already adjusted
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
        
        % save output data with units
        outputData.units = units;
%         save(outputDataFileList(f).name,'outputData','-v7.3');
        
        % update duration for splitting files
        durationAll = durationAll + outputData.duration;
    end
    
    disp('done replacing spike info')

%% combine all .nev files into one, merge all stim timings into 1 as well
    pwd = cd;
    cd(inputData.folderpath);

    nevFileList = dirSorted('*_spikes.nev');
    outputDataFileList = dirSorted('*_outputData.mat');
    
    NEV_all = [];
    stimInfo = [];
    
    for nf = 1:numel(nevFileList)
        disp(num2str(nf))
        NEV_single = openNEV('read',[inputData.folderpath nevFileList(nf).name],'nosave');
        load(outputDataFileList(nf).name);
        
        if(nf == 1)
            NEV_all = NEV_single;
            NEV_all.Filename = [nevFileList(1).name(1:underscoreIdx(2)),'_spikes_all_merged.nev'];
            
            stimInfo = outputData.stimInfo;
            stimInfo.chanSent = outputData.waveforms.chanSent;
            stimInfo.waveSent = outputData.waveforms.waveSent;
            stimInfo.parameters = outputData.waveforms.parameters;
            
            durationAll_points = NEV_single.MetaTags.DataDuration;
            durationAll_sec = NEV_single.MetaTags.DataDurationSec;
        else
            % merge meta info - DataDuration, DataDurationSec
            NEV_all.MetaTags.DataDuration = NEV_all.MetaTags.DataDuration + NEV_single.MetaTags.DataDuration;
            NEV_all.MetaTags.DataDurationSec = NEV_all.MetaTags.DataDurationSec + NEV_single.MetaTags.DataDurationSec;
                        
            % merge spike data -- update time stamp, combine all matrices
            NEV_single.Data.Spikes.TimeStamp = NEV_single.Data.Spikes.TimeStamp + durationAll_points;
            
            NEV_all.Data.Spikes.TimeStamp = [NEV_all.Data.Spikes.TimeStamp,NEV_single.Data.Spikes.TimeStamp];
            NEV_all.Data.Spikes.Electrode = [NEV_all.Data.Spikes.Electrode,NEV_single.Data.Spikes.Electrode];
            NEV_all.Data.Spikes.Unit = [NEV_all.Data.Spikes.Unit,NEV_single.Data.Spikes.Unit];
            NEV_all.Data.Spikes.Waveform = [NEV_all.Data.Spikes.Waveform,NEV_single.Data.Spikes.Waveform];
            
            % merge serial data -- update time stamps, combine all matrices
            NEV_single.Data.SerialDigitalIO.TimeStamp = NEV_single.Data.SerialDigitalIO.TimeStamp + durationAll_points;
            NEV_single.Data.SerialDigitalIO.TimeStampSec = NEV_single.Data.SerialDigitalIO.TimeStampSec + durationAll_sec;
            
            NEV_all.Data.SerialDigitalIO.TimeStamp = [NEV_all.Data.SerialDigitalIO.TimeStamp, NEV_single.Data.SerialDigitalIO.TimeStamp];
            NEV_all.Data.SerialDigitalIO.TimeStampSec = [NEV_all.Data.SerialDigitalIO.TimeStampSec, NEV_single.Data.SerialDigitalIO.TimeStampSec];
            NEV_all.Data.SerialDigitalIO.InsertionReason = [NEV_all.Data.SerialDigitalIO.InsertionReason, NEV_single.Data.SerialDigitalIO.InsertionReason];
            NEV_all.Data.SerialDigitalIO.UnparsedData = [NEV_all.Data.SerialDigitalIO.UnparsedData; NEV_single.Data.SerialDigitalIO.UnparsedData];
            
            % merge stim info from output data to separate merged struct
            stimInfo.stimOn = [stimInfo.stimOn; outputData.stimInfo.stimOn + durationAll_points];
            stimInfo.stimOff = [stimInfo.stimOff; outputData.stimInfo.stimOff + durationAll_points];
            stimInfo.chanSent = [stimInfo.chanSent; outputData.waveforms.chanSent];
            stimInfo.waveSent = [stimInfo.waveSent; outputData.waveforms.waveSent];
            
            % update duration all
            durationAll_points = durationAll_points + NEV_single.MetaTags.DataDuration;
            durationAll_sec = durationAll_sec + NEV_single.MetaTags.DataDurationSec;
        end
    end
    
    NEV_all.MetaTags.DataDuration = durationAll_points;
    NEV_all.MetaTags.DataDurationSec = durationAll_sec;
    
    % rewrite big nev file
    underscoreIdx = strfind(nevFileList(1).name,'_');
    saveNEV(NEV_all,[inputData.folderpath nevFileList(1).name(1:underscoreIdx(2)),'_spikes_all_merged.nev'],'noreport');
    % save merged stim info
    save([inputData.folderpath nevFileList(1).name(1:underscoreIdx(2)),'_spikes_all_merged_stimInfo.mat'],'stimInfo','-v7.3');
    
    disp('done merging into one .nev')
%% load that file into a cds -- flag to add all of the stim related data
    cds = commonDataStructure();
    flist = dir('*spikes_all_merged.nev');
    cds.file2cds([inputData.folderpath flist(1).name],inputData.task,inputData.ranBy,inputData.monkey,inputData.labnum,'recoverPreSync');

