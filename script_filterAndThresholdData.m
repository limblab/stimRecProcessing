%% this script filters and thresholds stimulation data
    clear
    pwd = cd;
    inputData.folderpath= 'C:\Users\jts3256\Desktop\Han_20180331_doublePulse\'; % must have \ at the end
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

%% load in *_merged-s, save the units to each nev file
    disp('started saving unit files')

    pwd = cd;
    cd(inputData.folderpath);
    outputDataFileList = dirSorted('*outputData.mat');
    nevFileList = dirSorted('*_merged-s.NEV*');
    NEVname = nevFileList(1).name; % grab the first one
    
    NEV_data = openNEV('read', [inputData.folderpath NEVname],'nosave');

    for f = 1:1%numel(outputDataFileList)
        % split back into individual files
        
        % undo any duration adding do to resets
    
        % load normal nev
    
        % replace spike data in the normal NEV
        
        % save normal nev with a new name
    end
    
    disp('done replacing spike info')

%% combine all .nev files into one
    

%% load that file into a cds


