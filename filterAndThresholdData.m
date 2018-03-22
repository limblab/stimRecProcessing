function [outputFigures, outputData ] = filterAndThresholdData(inputData)
    % function that loads a file and filters/thresholds the data
    
    outputData=[];
    
    % check to see if file exists
    if(~isfield(inputData,'filename') || ~isfield(inputData,'folderpath') || exist([inputData.folderpath,inputData.filename],'file')==0)
        error('file was not provided or does not exist');
    end
    
    % check to see if file is an ns5
    if(strcmpi(inputData.filename(end-2:end),'ns5') == 0)
        error('file is not an ns5');
    end
    
    % load map file
    if(~isfield(inputData,'mapFile'))
        error('no map file')
    end
    mapData=loadMapFile(inputData.mapFile(8:end));
    
    % merge map file labels with array name <- probably not necessary
    for i=1:size(mapData,1)
        mapData.label{i}=[inputData.array1(6:end),mapData.label{i}];
    end
    
    %establish trackers for information that is common to all files:
    eList={};
    posList=[];
    chList=[];
    
    % build filter
    if(~isfield(inputData,'filterParams'))
        [bAcausal,aAcausal] = butter(6,500/(30000/2),'high');
        filterParams.b = bAcausal;
        filterParams.a = aAcausal;
        filterParams.order = 6;
        filterParams.type = 'high';
    else
        filterParams = inputData.filterParams;
        bAcausal = inputData.filterParams.b;
        aAcausal = inputData.filterParams.a;
    end
    
    % get threshold multiplier
    if(~isfield(inputData,'thresholdMult'))
        error('no threshold multiplier provided');
    end
    thresholdMult = inputData.thresholdMult;
    
    % variables to store spike information   
    if(~isfield(inputData,'preOffset'))
        error('no pre spike offset provided')
    end
    preOffset = inputData.preOffset;
    if(~isfield(inputData,'postOffset'))
        error('no post spike offset provided');
    end
    postOffset = inputData.postOffset;
    if(~isfield(inputData,'maxAmplitude'))
        error('no maximum spike amplitude provided');
    end
    maxAmplitude = inputData.maxAmplitude; %uV
    lengthWave = preOffset+postOffset+1;
    numZeros = 200;
    
    % cds and extraseconds for merge purposes
    nevData = [];
    rawData = [];
    
    % set up artifact data
    if(~isfield(inputData,'artifactDataTime'))
        error('no time specified for the artifact data');
    end
    artifactDataTime = inputData.artifactDataTime; % in ms
    artifactData.t = zeros(1000,1);
    artifactDataIndex = 1;
     
    %% load file
    disp(['working on:', inputData.filename])
        
    NSx=openNSx('read', [inputData.folderpath,inputData.filename]);
    % remove needless spaces from the NSx.ElectrodesInfo.Label field
    for j = 1:numel(NSx.ElectrodesInfo)
        NSx.ElectrodesInfo(j).Label = strtrim(NSx.ElectrodesInfo(j).Label); % remove spaces
        NSx.ElectrodesInfo(j).Label(double(NSx.ElectrodesInfo(j).Label)==0) = []; % remove null values           
    end
    % get channels (electrode Id < 97)
    chanMask = [NSx.ElectrodesInfo.ElectrodeID]<97;
    % preset artifactData.artifact size based on num channels, subtract 1
    % for the time column
    artifactData.artifact = zeros(1000,sum(chanMask),artifactDataTime*2*30000/1000);
        
    %% find sync signal in analog data
    useSync=true;
    NSx_syncIdx=[];
    syncIdx=[];

    if(~isfield(inputData,'useSyncLabel') || isempty(inputData.useSyncLabel))
    %look for sync under the label sync
        for j=1:numel(NSx.ElectrodesInfo)
            syncIdx=find(strcmp(NSx.ElectrodesInfo(j).Label,'sync'));
            if ~isempty(syncIdx)
                NSx_syncIdx=j;
                syncName='sync';
            end
        end
        %if it wasn't called sync, try for matt's 'StimTrig' label:
        if isempty(syncIdx)
            for j=1:numel(NSx.ElectrodesInfo)
                syncIdx=find(strcmp(NSx.ElectrodesInfo(j).Label,'StimTrig'));
                if ~isempty(syncIdx)
                    NSx_syncIdx=j;
                    syncName='StimTrig';
                end
            end
        end
        %if we didn't find a sync channel, just look for ainp16
        if isempty(syncIdx)
            for j=1:numel(NSx.ElectrodesInfo)
                syncIdx=find(strcmp(NSx.ElectrodesInfo(j).Label,'ainp16'));
                if ~isempty(syncIdx)
                    useSync=false;
                    NSx_syncIdx=j;
                    syncName='ainp16';
                end
            end
        end

    else
        useSync=inputData.useSyncLabel;
        if useSync
            %find aIdx:
            for j=1:numel(NSx.ElectrodesInfo)
                syncIdx=find(strcmp(NSx.ElectrodesInfo(j).Label,'sync'));
                if ~isempty(syncIdx)
                    NSx_syncIdx=j;
                    syncName='sync';
                end
            end
        else
            for j=1:numel(NSx.ElectrodesInfo)
                syncIdx=find(strcmp(NSx.ElectrodesInfo(j).Label,'ainp16'));
                if ~isempty(syncIdx)
                    useSync=false;
                    NSx_syncIdx=j;
                    syncName='ainp16';
                end
            end
        end
    end
    
    %% if recover presync, append data, store where data was combined
    outputData.preSyncTimes = [];
    outputData.preSyncPoints = [];
    data = [];
    if(isfield(inputData,'recoverPreSync') && inputData.recoverPreSync)
        for NSx_idx = 1:numel(NSx.Data)
            data = [data,NSx.Data{NSx_idx}(:,:)];
        end
        
        NSx.Data = {};
        NSx.Data{1} = data;
        clear data
        
        outputData.preSyncTimes = NSx.MetaTags.DataDurationSec;
        outputData.preSyncPoints = NSx.MetaTags.DataPoints;
        
    else
        
    end
    
    
    %% use sync to get stim times:
    artifactDataPre.stimOn=find(diff(cdsTemp.analog{NSx_syncIdx}.(syncName)-mean(cdsTemp.analog{NSx_syncIdx}.(syncName))>3)>.5);
    stimOff=find(diff(cdsTemp.analog{NSx_syncIdx}.(syncName)-mean(cdsTemp.analog{NSx_syncIdx}.(syncName))<-3)>.5);
    artifactDataPre.stimOn = artifactDataPre.stimOn;
    artifactDataPre.stimOff=nan(size(artifactDataPre.stimOn));
    for j=1:numel(artifactDataPre.stimOn)
        if j<numel(artifactDataPre.stimOn)
            next=artifactDataPre.stimOn(j+1);
        else
            next=numel(cdsTemp.analog{NSx_syncIdx}.(syncName));
        end
        offIdx=stimOff(find((stimOff>artifactDataPre.stimOn(j)& stimOff<next),1,'first'));
        if ~isempty(offIdx)
            artifactDataPre.stimOff(j)=offIdx;
        end
    end
        
    %% fix stim times if more than one pulse sent per wave
    if(inputData.moreThanOnePulsePerWave)
        stimOnTemp = [];
        for stimIdx = 1:numel(artifactDataPre.stimOn)
            stimOnTemp = [stimOnTemp; artifactDataPre.stimOn(stimIdx) + 1/inputData.pulseFrequency*30000*(0:1:(inputData.numPulses-1))'];
        end
        artifactDataPre.stimOn = stimOnTemp;
        artifactDataPre.stimOff = artifactDataPre.stimOn + 2;
    end
    
        %% make waveforms sent file if it does not exist
        underscoreIdx = find(fileList(i).name=='_');
        [~,fname,~] = fileparts(fileList(i).name);
        waveformFilename = strcat(fileList(i).name(1:underscoreIdx(end)),'waveformsSent',fname(underscoreIdx(end):end),'.mat');
        if(exist(waveformFilename)~=0)
            load(waveformFilename)
        end
        
        if(~exist(waveformFilename)~=0 && ~(noSync || noSyncIntended))
            [~,fname,~] = fileparts(fileList(i).name);
            waveformFilename = strcat(fileList(i).name(1:underscoreIdx(end)),'waveformsSent',fname(underscoreIdx(end):end),'.mat');
            waveforms.chanSent = -1*ones(numel(artifactDataPre.stimOn),1);
            waveforms.waveSent = 1*ones(numel(artifactDataPre.stimOn),1);
            waveforms.parameters{1,1} = [];
            save(waveformFilename,'waveforms','-v7.3');
        end
        
        %% extract data from cdsTemp so that the rest of the code moves quicker - table operations are slow af
        cdsTempLFP = cdsTemp.lfp{:,:};
        % if the duke board is connected, then get the data from the correct analog pin
        if(inputData.dukeBoardChannel > 0)
            cdsTempLFP(:,inputData.dukeBoardChannel+1) = -1*cdsTemp.analog{1,1}.(inputData.dukeBoardLabel)*1000/100;
        end
        if(cdsTempLFP(end,1) > cdsTemp.meta.duration-0.01) % apparently need to remove data?
            timeRemove = cdsTempLFP(end,1) - cdsTemp.meta.duration + 0.2;
            pointsRemove = timeRemove*30000;
            cdsTempLFP = cdsTempLFP(1:end-pointsRemove,:);
        end
        cdsTempDuration = cdsTemp.meta.duration;
        
        maskKeep = ones(numel(artifactDataPre.stimOn),1);
        for counter = 1:numel(artifactDataPre.stimOn)
            if(artifactDataPre.stimOn(counter)/30000 > cdsTemp.meta.duration - 0.1)
                maskKeep(counter,1) = 0;
                % remove data during this 1ms to avoid filtering problems
                % later
                cdsTempLFP(artifactDataPre.stimOn(counter):artifactDataPre.stimOn(counter)+1*30,:) = 0;
            end
        end
        artifactDataPre.stimOn = artifactDataPre.stimOn(logical(maskKeep));
        
        %% fix stim times if joe sees an issue with the data file
        if(inputData.issueExists)
            % cable falling out, remove stim times that this corresponds
            % too
            maskArtifactKeep = ones(numel(artifactDataPre.stimOn),1);
            for artIdx = 1:numel(artifactDataPre.stimOn)
                artData = cdsTempLFP(artifactDataPre.stimOn(artIdx):artifactDataPre.stimOn(artIdx)+90,inputData.dukeBoardChannel+1);
                if(mean(artData) < -8000 || max(artData-mean(artData)) < 500)
                    maskArtifactKeep(artIdx) = 0;
                end
            end
            artifactDataPre.stimOn = artifactDataPre.stimOn(maskArtifactKeep==1);
            % fix waveforms sent file
            waveformFilename = strcat(fileList(i).name(1:underscoreIdx(end)),'waveformsSent',fname(underscoreIdx(end):end),'.mat');
            load(waveformFilename)
             
            waveforms.chanSent = waveforms.chanSent(maskArtifactKeep == 1);
            waveforms.waveSent = waveforms.waveSent(maskArtifactKeep == 1);
            save(waveformFilename,'waveforms','-v7.3');
        end
        
        %% merge cdsTemp with cds <- this one is everything but lfp + analog
        if(strcmp(inputData.task,'taskCObump'))
            if(i==1)
                % rewrite everything into cds -- which is really not a
                % commonDataStructure object but whatever
                cds.offsets = [0];
                cds.meta = cdsTemp.meta;
                cds.meta.hasLfp = 0;
                cds.meta.preOffset = preOffset;
                cds.meta.postOffset = postOffset;
                if(~isempty(cdsTemp.kin))
                    cds.kin = cdsTemp.kin;
                end
                cds.force = cdsTemp.force;
                cds.lfp = {};
                cds.emg = cdsTemp.emg;
                cds.analog = {};
                cds.stimOn = cdsTempLFP(artifactDataPre.stimOn,1); % use the lfp times because artifact data might start at non-zero value
                cds.timeOffset = cdsTempLFP(1,1);
                cds.stimOff = cdsTempLFP(artifactDataPre.stimOff,1);
%                 cds.stimOn = artifactDataPre.stimOn/30000;
%                 cds.stimOff = artifactDataPre.stimOff/30000;
                cds.triggers = cdsTemp.triggers;
                cds.units = cdsTemp.units;
                cds.trials = cdsTemp.trials;
                cds.aliasList = cdsTemp.aliasList;
                cds.operationLog = cdsTemp.operationLog;
                cds.meta.duration = 0;
            else     
                % port over kin data
                if(~isempty(cdsTemp.kin))
                    tempKin = cdsTemp.kin;
                    tempKin.t = tempKin.t + cds.meta.duration;
                    cdsKin = cds.kin;
                    cdsKin{end+1:end+size(tempKin,1),:} = tempKin{:,:};
                    cds.kin = cdsKin;
                end
                % port over stimOn and stimOff data
                cds.stimOn = [cds.stimOn;cdsTempLFP(artifactDataPre.stimOn,1) + cds.meta.duration];
                cds.stimOff = [cds.stimOff;cdsTempLFP(artifactDataPre.stimOff,1) + cds.meta.duration];
                % port over trial data
                if(~isempty(cdsTemp.trials))
                    tempTrials = cdsTemp.trials;
                    tempTrials.number = tempTrials.number + cds.trials.number(end);
                    tempTrials.startTime = tempTrials.startTime + cds.meta.duration;
                    tempTrials.endTime = tempTrials.endTime + cds.meta.duration;
                    tempTrials.tgtOnTime = tempTrials.tgtOnTime + cds.meta.duration;
                    tempTrials.goCueTime = tempTrials.goCueTime + cds.meta.duration;
                    tempTrials.bumpTime = tempTrials.bumpTime + cds.meta.duration;
                    tempTrials.stimTime = tempTrials.stimTime + cds.meta.duration;
                    cdsTrials = cds.trials;
                    cdsTrials(end+1:end+size(cdsTemp.trials,1),:) = tempTrials(:,:);
                    cds.trials = cdsTrials;
                end
                % update meta information
                cds.offsets(end+1,1) = cds.meta.duration;
                cds.meta.dataWindow(2) = cds.meta.dataWindow(2) + cdsTemp.meta.dataWindow(2);
                cds.meta.numTrials = cds.meta.numTrials + cdsTemp.meta.numTrials;
                cds.meta.numReward = cds.meta.numReward + cdsTemp.meta.numReward;
                cds.meta.numAbort = cds.meta.numAbort + cdsTemp.meta.numAbort;
                cds.meta.numFail = cds.meta.numFail + cds.meta.numFail;
                cds.meta.numIncomplete = cds.meta.numIncomplete + cds.meta.numIncomplete;
            end            
        else % ignore trial table completely
            if(i==1)
                % rewrite everything into cds -- which is really not a
                % commonDataStructure object but whatever
                cds.offsets = [0];
                cds.meta = cdsTemp.meta;
                cds.meta.hasLfp = 0;
                if(~isempty(cdsTemp.kin))
                    cds.kin = cdsTemp.kin;
                end
                cds.force = cdsTemp.force;
                cds.lfp = {};
                cds.emg = cdsTemp.emg;
                cds.analog = {};
                if(~noSync)
                    cds.stimOn = cdsTempLFP(artifactDataPre.stimOn,1);
                    cds.timeOffset = cdsTempLFP(1,1);
                    cds.stimOff = cdsTempLFP(artifactDataPre.stimOff,1);
                end
                cds.triggers = cdsTemp.triggers;
                cds.units = cdsTemp.units;
                cds.aliasList = cdsTemp.aliasList;
                cds.operationLog = cdsTemp.operationLog;
                cds.meta.duration = 0;
            else     
                % port over kin data
                if(~isempty(cdsTemp.kin))
                    tempKin = cdsTemp.kin;
                    tempKin.t = tempKin.t + cds.meta.duration;
                    cdsKin = cds.kin;
                    cdsKin{end+1:end+size(tempKin,1),:} = tempKin{:,:};
                    cds.kin = cdsKin;
                end
                % port over stimOn and stimOff data
                if(~noSync)
                    cds.stimOn = [cds.stimOn;cdsTempLFP(artifactDataPre.stimOn,1) + cds.meta.duration];
                    cds.stimOff = [cds.stimOff;cdsTempLFP(artifactDataPre.stimOff,1) + cds.meta.duration];  
                end
                % update meta information
                cds.offsets(end+1,1) = cds.meta.duration;
                cds.meta.dataWindow(2) = cds.meta.dataWindow(2) + cdsTemp.meta.dataWindow(2);
            end
        end
        
        clear cdsTemp

        % get template if applicable
        if(inputData.templateSubtract && exist(waveformFilename)~=0)
            % need to know stimulated channels
            numStimChan = numel(unique(waveforms.chanSent));
            numStims = zeros(numStimChan,1);
            for stimChanIdx = 1:numStimChan
                stimChans = unique(waveforms.chanSent);
                numStims(stimChanIdx) = sum(waveforms.chanSent == stimChans(stimChanIdx));
            end
            templateAll = zeros(numStimChan,96,templateSize);
            for ch = 2:size(cdsTempLFP,2)
                for stimuli = 1:numel(artifactDataPre.stimOn)+1
                    if(numel(artifactDataPre.stimOn)==0)
                        stimData = cdsTempLFP(:,ch);
                    elseif(stimuli == 1) % all data before first stim
                        stimData = cdsTempLFP(1:artifactDataPre.stimOn(stimuli),ch);
                    elseif(stimuli == numel(artifactDataPre.stimOn)+1) % all data after last stim artifact
                        stimData = cdsTempLFP(artifactDataPre.stimOn(stimuli-1)+5*30:end,ch);
                    else % data before ith stim up to ith-1 
                        if(artifactDataPre.stimOn(stimuli) - 6*30 > artifactDataPre.stimOn(stimuli-1))  
                            try
                                stimData = cdsTempLFP(artifactDataPre.stimOn(stimuli-1):artifactDataPre.stimOn(stimuli),ch);
                                stimChan = find(unique(waveforms.chanSent) == waveforms.chanSent(stimuli-1));
                                stimDataTemp = [stimData(:,1);zeros(numZeros,1)];
                                stimDataTemp = fliplr(filter(bAcausal,aAcausal,fliplr(stimDataTemp')))';
                                stimDataTemp = stimDataTemp(1:end-numZeros,1);
                                templateAll(stimChan,ch-1,:) = squeeze(templateAll(stimChan,ch-1,:)) + (stimDataTemp(1:templateSize,1)/numStims(stimChan));
                            catch
                                disp('here');
                            end
                            
                        end
                    end
                    
                end
            end
            
        end
        
        %% add fake stim times so that the data processed is not too large and does not slow things down
        artifactDataPreMask = ones(numel(artifactDataPre.stimOn),1);
        artifactDataPreStimOnTemp = artifactDataPre.stimOn;
        tempIdx = 2;
        while tempIdx < numel(artifactDataPreStimOnTemp)
            if(artifactDataPreStimOnTemp(tempIdx) - artifactDataPreStimOnTemp(tempIdx-1) > inputData.maxChunkLength)
                artifactDataPreStimOnTemp = [artifactDataPreStimOnTemp(1:tempIdx-1,1);artifactDataPreStimOnTemp(tempIdx-1)+inputData.maxChunkLength;...
                    artifactDataPreStimOnTemp(tempIdx:end,1)];
                artifactDataPreMask = [artifactDataPreMask(1:tempIdx-1,1); 0; artifactDataPreMask(tempIdx:end,1)];
            end
            tempIdx = tempIdx + 1;
        end
        artifactDataPre.stimOn = artifactDataPreStimOnTemp;
        
        
        %% get thresholds for each channel based on non stim data
        thresholdAll = zeros(size(cdsTempLFP,2)-1,1);
        flagStop = 0;
        for ch = 2:size(cdsTempLFP,2)
            numPoints = 0;
            for stimuli = 1:numel(artifactDataPre.stimOn)+1
                if(numel(artifactDataPre.stimOn)==0)
                    stimData = cdsTempLFP(:,ch);
                elseif(stimuli == 1) % all data before first stim
                    stimData = cdsTempLFP(1:artifactDataPre.stimOn(stimuli)-1*30,ch);
                elseif(stimuli == numel(artifactDataPre.stimOn) + 1) % all data after last stim artifact
                    stimData = cdsTempLFP(artifactDataPre.stimOn(stimuli-1)+5*30:end,ch);
                else % data before ith stim up to ith-1 
                    if(artifactDataPre.stimOn(stimuli) - 5*30 > artifactDataPre.stimOn(stimuli-1))  
                        stimData = cdsTempLFP(artifactDataPre.stimOn(stimuli-1)+5*30:artifactDataPre.stimOn(stimuli)-1*30,ch);
                    end
                end
                try
                    stimDataTemp = [stimData(:,1);mean(stimData(end-20:end,1))*ones(numZeros,1)];
                    numPoints = numPoints + numel(stimData);
                    stimDataTemp = fliplr(filter(bAcausal,aAcausal,fliplr(stimDataTemp')))';
                    stimDataTemp = stimDataTemp(1:end-numZeros,1);
                    thresholdAll(ch-1) = thresholdAll(ch-1) + sum(stimDataTemp.^2);
                catch
                end
            end
            thresholdAll(ch-1) = sqrt(thresholdAll(ch-1)/numPoints);
        end
        
        
        spikeWaves = zeros(10000,lengthWave);
        spikeTimes = zeros(10000,1);
        spikeChan = zeros(10000,1);    
        spikeNum = 1;
        disp('filtering and thresholding')
        for ch = 1:(size(cdsTempLFP,2)-1)
            for stimIdx = 1:numel(artifactDataPre.stimOn)+1
                stimData = [];
                if(numel(artifactDataPre.stimOn) == 0)
                    stimData = cdsTempLFP(:,2:end);
                elseif(stimIdx == 1) % all data before first stim
                    stimData = cdsTempLFP(1:artifactDataPre.stimOn(stimIdx),2:end);
                elseif(stimIdx == numel(artifactDataPre.stimOn) + 1) % all data after last stim
                    stimData = cdsTempLFP(artifactDataPre.stimOn(stimIdx-1):end,2:end);
    %                 % perform pca step
    %                 [coeff,score] = pca(stimData);
    %                 coeff(:,1:4) = 0;
    %                 stimData = coeff*score';
    %                 stimData = stimData';
                    stimData(1:max(1,inputData.blankPeriod),:) = 0; % blank the first millisecond
                else % data before ith stim up to ith-1 
                    stimData = cdsTempLFP(artifactDataPre.stimOn(stimIdx-1):artifactDataPre.stimOn(stimIdx),2:end);
                    % perform pca step
    %                 [coeff,score] = pca(stimData);
    %                 coeff(:,1:4) = 0;
    %                 stimData = coeff*score';
    %                 stimData = stimData';
                    stimData(1:max(1,inputData.blankPeriod),:) = 0; % blank 
                end
            
                % filter backwards on all channels and threshold
                stimDataTemp = [stimData(:,ch);mean(stimData(end-min(length(stimData(:,ch))-1,20):end,ch))*ones(numZeros,1)];
                
                stimDataTemp = fliplr(filter(bAcausal,aAcausal,fliplr(stimDataTemp')))';
                stimData(:,ch) = stimDataTemp(1:end-numZeros);
                if(inputData.templateSubtract && stimIdx~=1)
                    stimChan = find(unique(waveforms.chanSent) == waveforms.chanSent(stimIdx-1));
                    if(size(stimData,1) > size(templateAll,3))
                        stimData(1:templateSize,ch) = stimData(1:templateSize,ch) - squeeze(templateAll(stimChan,ch,:));
                    end
                end
                if(inputData.templateSubtract && any(isfield(inputData,'blankPeriod')))
                    stimData(1:max(1,inputData.blankPeriod),ch) = 0;
                else
%                     stimData(1:30*1,ch) = 0;
                end
                threshold = thresholdMult*thresholdAll(ch);
                if(abs(threshold) < 1)
                    threshold = sign(threshold)*100000;
                end
                
                %% get threshold crossings
%                 threshold = abs(rms(stimData(max(1,numel(stimData(:,ch))-10):end,ch))*thresholdMult);
                thresholdCrossings = find(stimData(:,ch)>abs(threshold));
                
                %% append data before and after stimData to get spikes near the edges
                numAppend = 100;
                stimData = [zeros(numAppend,size(stimData,2));stimData(:,:);zeros(numAppend,size(stimData,2))];
                thresholdCrossings = thresholdCrossings + numAppend;
                
                % remove potential artifacts and too close to beginning/end
                % of stim data
                crossingsMask = ones(numel(thresholdCrossings),1);
                for cross = 1:numel(thresholdCrossings)
                    if(stimData(thresholdCrossings(cross),ch) > maxAmplitude || ...
                            ~(thresholdCrossings(cross)+postOffset <= numel(stimData(:,ch)) && thresholdCrossings(cross)-preOffset > 0))
                        crossingsMask(cross) = 0;
                    end
                end
                thresholdCrossings = thresholdCrossings(crossingsMask(:) == 1);
                
                
                % remove chains -- find largest spot
                idx = 2;
                chain = [1];
                crossingsKeep = [];
                while idx <= numel(thresholdCrossings)
                    if(thresholdCrossings(idx) == thresholdCrossings(idx-1)+1) % store in chain
                        chain = [chain;idx];
                    elseif(~isempty(chain)) % broke a chain, store minidx, update idx, empty chain
                        [~,maxIdx] = max(stimData(thresholdCrossings(chain)));
                        if(isempty(crossingsKeep))
                            crossingsKeep = [thresholdCrossings(maxIdx+chain(1)-1)];
                        else
                            crossingsKeep = [crossingsKeep;thresholdCrossings(maxIdx+chain(1)-1)];
                        end
                        chain = [idx];
                    end
                    idx = idx+1;
                end
                if(numel(thresholdCrossings) > 0)
                    thresholdCrossings = [crossingsKeep;thresholdCrossings(end)];
                end
                
                % go through and weed out ones that are too close to each other
                % prioritize backwards in time
                crossingsMask = ones(numel(thresholdCrossings),1);
                for cross = numel(thresholdCrossings):-1:2
                    if(crossingsMask(cross) == 1) % check time beforehand to see if one is too close
                        crossCheck = cross-1;
                        while crossCheck >= 1 && thresholdCrossings(crossCheck) >= thresholdCrossings(cross) - max(preOffset,postOffset)
                            crossingsMask(crossCheck) = 0;
                            crossCheck = crossCheck-1;
                        end
                    end
                end  
                thresholdCrossings = thresholdCrossings(crossingsMask(:) == 1);

                % store thresholdCrossing data
                for cross = 1:numel(thresholdCrossings)
                    if(stimIdx == 1)
                        spikeTimes(spikeNum) = cdsTempLFP(thresholdCrossings(cross)-numAppend,1); % this is in seconds
                        
                    else
                        spikeTimes(spikeNum) = cdsTempLFP(artifactDataPre.stimOn(stimIdx-1)+thresholdCrossings(cross)-1-numAppend,1); % this is in secondss
                    end
                    spikeChan(spikeNum) = ch;
                    % check if too close to beginning or end
                    if(thresholdCrossings(cross)+postOffset > numel(stimData(:,ch)))
                        numZerosPad = (preOffset+postOffset+1) - numel(stimData(thresholdCrossings(cross)-preOffset:end,ch));
                        spikeWaves(spikeNum,:) = [stimData(thresholdCrossings(cross)-preOffset:end,ch);zeros(numZerosPad,1)];
                    elseif(thresholdCrossings(cross)-preOffset < 0)
                        numZerosPad = (preOffset+postOffset+1) - numel(stimData(1:thresholdCrossings(cross)+postOffset,ch));
                        spikeWaves(spikeNum,:) = [zeros(numZerosPad,1);stimData(thresholdCrossings(cross)-preOffset:end,ch)];
                    else
                        spikeWaves(spikeNum,:) = stimData(thresholdCrossings(cross)-preOffset:thresholdCrossings(cross)+postOffset,ch);
                    end
                    spikeNum = spikeNum + 1;

                    if(spikeNum > 0.67*numel(spikeTimes))
                        spikeTimes = [spikeTimes;zeros(1000,1)];
                        spikeWaves = [spikeWaves; zeros(1000,lengthWave)];
                        spikeChan = [spikeChan; zeros(1000,1)];
                    end
                end
                 
                
            end % end channel for
        end % end stimOn for
        
        % store spike data
        disp('storing data')
        spikeTimes = spikeTimes(1:spikeNum-1,1);
        spikeWaves = spikeWaves(1:spikeNum-1,:);
        spikeChan = spikeChan(1:spikeNum-1,:);
        [spikeTimes,sortOrder] = sort(spikeTimes);
        spikeWaves = spikeWaves(sortOrder,:);
        spikeChan = spikeChan(sortOrder,:);
        if(i==1)
            nevData.ts = spikeTimes;
            nevData.waveforms = spikeWaves(:,:);
            nevData.elec = spikeChan(:,:);
        else
            nevData.ts(end+1:end+numel(spikeTimes),:) = spikeTimes + cds.meta.duration;
            nevData.waveforms(end+1:end+numel(spikeTimes),:) = spikeWaves(:,:);
            nevData.elec(end+1:end+numel(spikeTimes),:) = spikeChan(:,:);
        end
        % store raw data where a threshold crossing occurs (index of
        % spikeTimes in this file)
        if(i==1)
            rawData.ts = spikeTimes;
            rawData.waveforms = zeros(numel(spikeTimes),preOffset*2+postOffset*2+1);
            rawIdxList = zeros(numel(spikeTimes),1);
            stIdx = 1;
            for ttt = 1:size(cdsTempLFP,1)
                while(stIdx <= numel(spikeTimes) && cdsTempLFP(ttt,1) == spikeTimes(stIdx))
                    rawIdxList(stIdx,1) = ttt;
                    stIdx = stIdx + 1;
                end
            end
            for r = 1:numel(spikeTimes)
                rawIdx = rawIdxList(r);
                if(rawIdx ~= 0 && rawIdx-preOffset*2 > 0 && rawIdx+postOffset*2 <= size(cdsTempLFP,1))
                    rawData.waveforms(r,:) = cdsTempLFP(rawIdx-preOffset*2:rawIdx+postOffset*2,spikeChan(r)+1);
                end
            end
            rawData.elec = spikeChan;
        else
            rawData.ts(end+1:end+numel(spikeTimes),:) = spikeTimes + cds.meta.duration;
            rOffset = size(rawData.waveforms,1);
            rawData.waveforms = zeros(numel(spikeTimes),preOffset*2+postOffset*2+1);
            rawIdxList = zeros(numel(spikeTimes),1);
            stIdx = 1;
            for ttt = 1:size(cdsTempLFP,1)
                while(stIdx <= numel(spikeTimes) && cdsTempLFP(ttt,1) == spikeTimes(stIdx))
                    rawIdxList(stIdx,1) = ttt;
                    stIdx = stIdx + 1;
                end
            end
            for r = 1:numel(spikeTimes)
                rawIdx = rawIdxList(r);
                if(rawIdx ~= 0 && rawIdx-preOffset*2 > 0 && rawIdx+postOffset*2 <= size(cdsTempLFP,1))
                    rawData.waveforms(r+rOffset,:) = cdsTempLFP(rawIdx-preOffset*2:rawIdx+postOffset*2,spikeChan(r)+1);
                end
            end
            rawData.elec = spikeChan;
        end

        clear spikeTimes
        clear spikeChan
        clear spikeWaves
        
        % store artifact data
        artifactDataPre.stimOn = artifactDataPre.stimOn(artifactDataPreMask==1,1);
        for art = 1:inputData.artifactSkip:numel(artifactDataPre.stimOn)           
            if(artifactDataPre.stimOn(art) + artifactDataTime*30000/1000 <= size(cdsTempLFP,1))
                artifactData.artifact(artifactDataIndex,:,:) = cdsTempLFP(artifactDataPre.stimOn(art)-floor(artifactDataTime*30000/1000):artifactDataPre.stimOn(art)+floor(artifactDataTime*30000/1000)-1,2:end)';
                if(i==1)
                    artifactData.t(artifactDataIndex,1) = cdsTempLFP(artifactDataPre.stimOn(art),1);
                else
                    artifactData.t(artifactDataIndex,1) = cdsTempLFP(artifactDataPre.stimOn(art),1) + cds.meta.duration;
                end
                artifactDataIndex = artifactDataIndex + 1;
            end
            if(artifactDataIndex >= size(artifactData.artifact,1))
                artifactData.artifact(end+1:end+1000,:,:) = 0;
                artifactData.t(end+1:end+1000,1) = 0;
            end
        end
        
%         cds.meta.duration = cds.meta.duration + cdsTempDuration;
       
%     cds.units = [];
%     cds.artifactData = artifactData;
%     cds.rawData = rawData;
%     if(any(isfield(inputData,'stimsPerBump')))
%         cds.stimsPerBump = inputData.stimsPerBump;
%     else
%         cds.stimsPerBump = -1; % usually means no bump, can also mean that the user forgot to input stimsPerBump
%     end
%     cds.filter = filterParams;
%     cds.processed = 0;
%     save(strcat(fileList(1).name(1:end-4),'_cds.mat'),'cds','-v7.3');
%     save(strcat(fileList(1).name(1:end-4),'_nevData.mat'),'nevData','-v7.3');
end