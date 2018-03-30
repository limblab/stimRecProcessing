%% this script filters and thresholds stimulation data
    clear
    pwd = cd;
    inputData.folderpath= 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\testingCode\'; % must have \ at the end
    inputData.mapFile='mapFileR:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
    % inputData.mapFile = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Chips_12H1\map_files\left S1\SN 6251-001455.cmp';

    inputData.task='taskCO';
    inputData.ranBy='ranByJoseph'; 
    inputData.array1='arrayLeftS1'; 
    inputData.monkey='monkeyChips';
    iniputData.labnum = 6;
    
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
    filename = strcat(fileList(1).name(1:end-4),'_merged');
    mapFilename = inputData.mapFile(8:end);
    comments = '';
    writeNEV(nevDataAll, packetWidth, filename, mapFilename, comments )

    cd(pwd);
    disp('DONE -- CAN CONTINUE')
    warning('on')

%% sort *_merged and call it *_merged-s

%% load in *_merged-s, save the units to each file
    disp('started saving unit files')

    pwd = cd;
    cd(inputData.folderpath);
    outputDataFileList = dirSorted('*outputData.mat');
    nevFileList = dirSorted('*_merged.NEV*');
    nevFile = nevFileList(1).name; % grab the first one
    

    spikes = commonDataStructure();
    spikes.file2cds([inputData.folderpath,nevFileList(1).name],inputData.ranBy,inputData.array1,...
        inputData.monkey,inputData.labnum,'ignoreJumps',inputData.task,inputData.mapFile);
    
    % split spikes units and store in a separate file
    durationAll = 0;
    for f = 1:numel(outputDataFileList)
        load(outputDataFileList(f).name);
        units = spikes.units; % move out of a cds so it can be changed
        
        for nn = 1:size(spikes.units,2)
            units(nn).spikes.ts = units(nn).spikes.ts-durationAll;
            spikeMask = units(nn).spikes.ts > 0 & ...
                units(nn).spikes.ts < outputData.duration;
            units(nn).spikes = units(nn).spikes(spikeMask,:);
        end
        outputData.units = units;
        save(strcat(outputDataFileList(f).name(1:end-4),'_units'),'outputData');
        durationAll = durationAll + outputData.duration;
    end
    disp('done with this step')

%% for each .ns5 and .nev, merge unit data from above with all the other data into a single .nev file
    disp('started merging unit data with nev data')
    pwd = cd;
    cd(inputData.folderpath);
    outputDataFileList = dirSorted('*_units.mat');
    
    for f = 1:numel(outputDataFileList)
        % load in output data
        load(outputDataFileList(f).name);
        % load in all corresponding nev files
        NEV_data = nev2NEVNSx_noNS5(outputDataFileList(f).name(1:end-3));
        
        % merge magically
        
    end


%% combine all .nev files
    



% %% merge spike and stimulation data into one file with waveforms sent information
% % folderpath = 'R:\data\Mihili_12A3\stimRecord\Mihili_20170717_stimRecord\';
% pwd = cd;
% cd(folderpath)
% fileListProcessed = dirSorted('*processed.mat');
% dataAll = [];
% for f = 1:numel(fileListProcessed)
%     load(fileListProcessed(f).name);
%     
%     % merge units, waveforms sent, raw data, and artifact data into one struct,
%     % ditch nev data
%     
%     if(f==1)
%         % copy everything over
%         dataAll.units = outputData.units;
%         dataAll.waveforms = outputData.waveforms;
%         dataAll.waveforms.parameters{1}.numStims = numel(outputData.waveforms.chanSent);
%         dataAll.rawData = outputData.rawData;
%         dataAll.artifactData = outputData.artifactData;
%         dataAll.duration = outputData.duration;
%     else     
%         % update units information - i think this is covering for the case
%         % when not all units show up in a file? but poorly
%         for nn = 1:size(outputData.units,2)
%             nnAll = 1;
%             outputData.units(nn).spikes{:,1} = outputData.units(nn).spikes{:,1} + dataAll.duration;
%             while(nnAll <= size(outputData.units,2) && (outputData.units(nnAll).chan ~= outputData.units(nn).chan || outputData.units(nnAll).ID ~= outputData.units(nn).ID))
%                 nnAll = nnAll + 1;
%             end
%             if(nnAll <= size(cds.units,2))
%                 dataAll.units(nnAll).spikes{end+1:end+size(cds.units(nn).spikes,1),:} = outputData.units(nn).spikes{:,:};
%             end
%             
%         end
%         
%         % update waveforms data
%         dataAll.waveforms.chanSent(end+1:end+numel(outputData.waveforms.chanSent)) = outputData.waveforms.chanSent;
%         dataAll.waveforms.waveSent(end+1:end+numel(outputData.waveforms.waveSent)) = outputData.waveforms.waveSent;
%         dataAll.waveforms.parameters{end+1} = outputData.waveforms.parameters;
%         dataAll.waveforms.parameters{end+1}.numStims = numel(outputData.waveforms.chanSent);
% 
%         % update raw data
%         dataAll.rawData.ts(end+1:end+numel(outputData.rawData.ts)) = outputData.rawData.ts;
%         dataAll.rawData.ts(end+1:end+numel(outputData.rawData.ts),:) = outputData.rawData.waveforms;
%         dataAll.rawData.ts(end+1:end+numel(outputData.rawData.elec)) = outputData.rawData.elec;
%         
%         % update artifact data
%         dataAll.artifactData.t(end+1:end+numel(outputData.artifactData.t)) = outputData.artifactData.t;
%         dataAll.artifactData.artifact(end+1:end+numel(outputData.artifactData.t),:,:) = outputData.artifactData.artifact;
%         
%         % update duration 
%         dataAll.duration = dataAll.duration + outputData.duration;
%     end
%     
% end
% 
% outputData = dataAll;
% clear dataAll
% save(strcat(fileListProcessed(1).name(1:26),'_all_processed'),'outputData','-v7.3');
% 
% cd(pwd)
% disp('done merging')
