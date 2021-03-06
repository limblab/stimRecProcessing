function [outputFigures, outputData ] = filterAndThresholdData(inputData)
    % function that loads a file and filters/thresholds the data
    
    outputData=[];
    outputFigures = [];
    

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
    
    % check to see if good chan list exists
    if(~isfield(inputData,'good_chan_list'))
        error('need to supply good_chan_list')
    end
    
    mapData=loadMapFile(inputData.mapFile(8:end));
    
    % check if presample exists
    if(~isfield(inputData,'presample'))
        inputData.presample = 100;
        warning('no presample provided, defaulting to 100');
    end
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
    
    % initialize arrays
    nevData = [];
    rawData = [];
    
    % check for artifact data time
    if(~isfield(inputData,'artifactDataTime'))
        error('no time specified for the artifact data');
    end
    artifactDataTime = inputData.artifactDataTime; % in ms
     
    if(~isfield(inputData,'softwareAmplification'))
        error('no softwareAmplification field');
    elseif(inputData.softwareAmplification)
        % try to find file with glm_data in filename
        gain_ratio_file = dir([inputData.folderpath,'*gain_ratio*']);
        if(~isempty(gain_ratio_file))
            load([gain_ratio_file.folder,'/',gain_ratio_file.name]);
            inputData.dukeBoard.gain_ratio = gain_ratio;
            inputData.dukeBoard.amp_1_all = amp_1_all;
            inputData.dukeBoard.polarity = polarity;
            inputData.dukeBoard.t_post_stim = t_post_stim_test;
            inputData.dukeBoard.pw1 = pw1;
            clear gain_ratio; clear amp_1_all; clear polarity; clear t_post_stim_test;
            clear pw1;
            
            warning('loaded a gain_ratio from a file. Stop if this is not desired');
        else
            % else throw an error
            error('supposed to do software amplification but glm field is empty and can''t find file with glm_data in name');
        end
    end
        
    %% load file
    disp(['working on:', inputData.filename])

    try
        NSx=openNSx('read', [inputData.folderpath,inputData.filename],'precision','double','uV');
    catch
        NSx=openNSx('read', [inputData.folderpath,inputData.filename],'precision','double');
        warning('openNSx does not support the uV option (?). This may somehow affect the amplitude of collected spikes');
    end
    
    NSx_trim = []; % store a version for trimming
    
    
    %% remove needless spaces from the NSx.ElectrodesInfo.Label field

%     for j = 1:numel(NSx.ElectrodesInfo)
%         NSx.ElectrodesInfo(j).Label = strtrim(NSx.ElectrodesInfo(j).Label); % remove spaces
%         NSx.ElectrodesInfo(j).Label(double(NSx.ElectrodesInfo(j).Label)==0) = []; % remove null values           
%     end
    % get channels (electrode Id < 97)
    chanMask = [NSx.ElectrodesInfo.ElectrodeID]<97;
    chanMapping = [NSx.ElectrodesInfo.ElectrodeID];    

    % map duke board channel to correct electrode, remove electrode from
    % chanMask
    if(isfield(inputData,'dukeBoardChannel') && isfield(inputData,'dukeBoardLabel') && inputData.dukeBoardChannel > 0)
        for j = 1:numel(NSx.ElectrodesInfo)
            if(~isempty(strfind(NSx.ElectrodesInfo(j).Label,inputData.dukeBoardLabel)))
                % make chanMask and chanMapping correct
                chanMask(j) = 1;
                chanMapping(j) = inputData.dukeBoardChannel;
                % remove data recorded on the electrode
                elecIdx = find([NSx.ElectrodesInfo.ElectrodeID] == inputData.dukeBoardChannel);
                chanMask(elecIdx) = 0;
            end
        end
    end
    
    %% find sync signal in analog data
    useSync=true;
    NSx_syncIdx=[];
    syncIdx=[];

    if(~isfield(inputData,'useSyncLabel') || isempty(inputData.useSyncLabel))
    %look for sync under the label sync
        for j=1:numel(NSx.ElectrodesInfo)
            syncIdx=strfind(NSx.ElectrodesInfo(j).Label,'sync');
            if ~isempty(syncIdx)
                NSx_syncIdx=j;
                syncName='sync';
            end
        end
        %if it wasn't called sync, try for matt's 'StimTrig' label:
        if isempty(syncIdx)
            for j=1:numel(NSx.ElectrodesInfo)
                syncIdx=strfind(NSx.ElectrodesInfo(j).Label,'StimTrig');
                if ~isempty(syncIdx)
                    NSx_syncIdx=j;
                    syncName='StimTrig';
                end
            end
        end
        %if we didn't find a sync channel, just look for ainp16
        if isempty(syncIdx)
            for j=1:numel(NSx.ElectrodesInfo)
                syncIdx=strfind(NSx.ElectrodesInfo(j).Label,'ainp16');
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
                syncIdx=strfind(NSx.ElectrodesInfo(j).Label,'sync');
                if ~isempty(syncIdx)
                    NSx_syncIdx=j;
                    syncName='sync';
                end
            end
        else
            for j=1:numel(NSx.ElectrodesInfo)
                syncIdx=strfind(NSx.ElectrodesInfo(j).Label,'ainp16');
                if ~isempty(syncIdx)
                    useSync=false;
                    NSx_syncIdx=j;
                    syncName='ainp16';
                end
            end
        end
    end
    
    %% remove channels from ns5, write back as a new ns5 for later use
    % this is so that the time stamps of non-neural data can be adjusted in the same way
    % as all of the other data
    fieldNames = fieldnames(NSx);
    for field_idx = 1:numel(fieldNames)
        if(strcmp(fieldNames(field_idx),'Data') == 1)
            if(iscell(NSx.Data))
                for data_idx = 1:numel(NSx.Data)
                    NSx_trim.Data{data_idx} = NSx.Data{data_idx}(~chanMask,:);
                end
            else
                NSx_trim.Data = NSx.Data(~chanMask,:);
            end
        else
            NSx_trim.(fieldNames{field_idx}) = NSx.(fieldNames{field_idx});
        end
    end
    NSx_trim.ElectrodesInfo(chanMask) = [];
    
    saveNSx(NSx_trim,[inputData.folderpath,inputData.filename(1:end-4) '_spikesExtracted.ns5'],'noreport');
   
    %% append data, store where data was combined
    
    outputData.preSyncTimes = [];
    outputData.preSyncPoints = [];
    data = [];
    if(iscell(NSx.Data))
        for NSx_idx = 1:numel(NSx.Data)
            if(NSx_idx == 1)
                data = NSx.Data{NSx_idx};
            else
                % to line up with the cds
                data = [data,zeros(size(NSx.Data{NSx_idx},1),NSx.MetaTags.Timestamp(NSx_idx)),NSx.Data{NSx_idx}(:,:)];
            end
        end 
        NSx.Data = data; 
        clear data
        
    end
        
    outputData.DataDurationSec = NSx.MetaTags.DataDurationSec + NSx.MetaTags.Timestamp/30000;
    outputData.DataPoints = NSx.MetaTags.DataPoints + NSx.MetaTags.Timestamp;
    outputData.TimeStamp = NSx.MetaTags.Timestamp;
    
    NSx_dataIdx = 1;
    outputData.duration = size(NSx.Data,2)/30000; % 

    
    %% use sync to get stim times:
    stimulationInformation.stimOn=find(diff(NSx.Data(NSx_syncIdx,:)-mean(NSx.Data(NSx_syncIdx,:))>3)>.5);
    stimOff=find(diff(NSx.Data(NSx_syncIdx,:)-mean(NSx.Data(NSx_syncIdx,:))<-3)>.5);
    stimulationInformation.stimOff=nan(size(stimulationInformation.stimOn));
    for j=1:numel(stimulationInformation.stimOn)
        if j<numel(stimulationInformation.stimOn)
            next=stimulationInformation.stimOn(j+1);
        else
            next=numel(NSx.Data(NSx_syncIdx,:));
        end
        offIdx=stimOff(find((stimOff>stimulationInformation.stimOn(j)& stimOff<next),1,'first'));
        if ~isempty(offIdx)
            stimulationInformation.stimOff(j)=offIdx;
        end
    end
    
    stimulationInformation.stimOn = stimulationInformation.stimOn';
    stimulationInformation.stimOff = stimulationInformation.stimOff';
        
    %% fix stim times if more than one pulse sent per wave
    if(isfield(inputData,'moreThanOnePulsePerWave') && isfield(inputData,'pulseFrequency') && isfield(inputData,'numPulses') && inputData.moreThanOnePulsePerWave)
        stimOnTemp = [];
        for stimIdx = 1:numel(stimulationInformation.stimOn)
            stimOnTemp = [stimOnTemp; stimulationInformation.stimOn(stimIdx) + floor(1/inputData.pulseFrequency*30000*(0:1:(inputData.numPulses-1)))'];
        end
        stimulationInformation.stimOn = stimOnTemp;
        stimulationInformation.stimOff = stimulationInformation.stimOn + 2;
    else
        warning('not handling more than one pulse per sync pulse because all inputs were not provided or this was flagged false');
    end
    
    %% make waveforms sent file if it does not exist
    underscoreIdx = find(inputData.filename=='_');
    [~,fname,~] = fileparts(inputData.filename);
    waveformFilename = strcat(fname(1:underscoreIdx(end)),'waveformsSent',fname(underscoreIdx(end):end),'.mat');
    waveformFilenameAbridged = strcat(fname(1:underscoreIdx(end)),'waveformsSent_',num2str(str2num(fname(underscoreIdx(end)+1:end))),'.mat');
    if(exist([inputData.folderpath,waveformFilename])~=0)
        load([inputData.folderpath,waveformFilename])
    elseif(exist([inputData.folderpath,waveformFilenameAbridged])~=0)
        load([inputData.folderpath,waveformFilenameAbridged]);
    elseif(~exist(waveformFilename)~=0)
        warning('no waveform sent file exists, continuing with an arbitrary one');
        waveformFilename = strcat(fname(1:underscoreIdx(end)),'waveformsSent',fname(underscoreIdx(end):end),'.mat');
        waveforms.chanSent = -1*ones(numel(stimulationInformation.stimOn),1);
        waveforms.waveSent = 1*ones(numel(stimulationInformation.stimOn),1);
        waveforms.parameters{1,1} = [];
        save([inputData.folderpath,waveformFilename],'waveforms','-v7.3');
    end
        
    %% add fake stim times so that the data processed is not too large and does not slow things down
    %% this might not need to be here anymore because ... (oops) 
    if(isfield(inputData,'maxChunkLength'))
        stimOnMask = ones(numel(stimulationInformation.stimOn),1);
        stimOnTemp = stimulationInformation.stimOn;
        tempIdx = 2;
        while tempIdx < numel(stimOnTemp)
            if(stimOnTemp(tempIdx) - stimOnTemp(tempIdx-1) > inputData.maxChunkLength + 1000*30) % 1 second extra
                stimOnTemp = [stimOnTemp(1:tempIdx-1,1);stimOnTemp(tempIdx-1)+inputData.maxChunkLength;...
                    stimOnTemp(tempIdx:end,1)];
                stimOnMask = [stimOnMask(1:tempIdx-1,1); 0; stimOnMask(tempIdx:end,1)];
            end
            tempIdx = tempIdx + 1;
        end
        stimulationInformation.stimOn = stimOnTemp;
    else
        warning('no max processing size provided, function may run slowly');
    end
        
        
    %% extract neural channels from NSx.Data, convert to double, and append time information
    %% this actually takes a decent chunk of time, room for improvement probably
    neuralLFP = zeros(sum(chanMask)+1,size(NSx.Data,2)); % preallocate
    
    neuralLFP(1,:) = roundTime((0:size(neuralLFP,2)-1)/NSx.MetaTags.SamplingFreq) + NSx.MetaTags.Timestamp(NSx_dataIdx)/NSx.MetaTags.TimeRes; % time stamps
    for ch = 1:size(NSx.Data,1)
        if(chanMask(ch))
            neuralLFP(chanMapping(ch)+1,:) = double(NSx.Data(ch,:)); % idx 1 is for time
        end
    end
    % NSx.Data = [];
    
    %% get thresholds for each channel based on non stim data. 
    thresholdAll = zeros(size(neuralLFP,1)-1,1);
    
    numPoints = size(neuralLFP,2);
    disp('thresholding data')
    numPointsActual = 0;
    for stimuli = 1:numel(stimulationInformation.stimOn)+1
        if(numel(stimulationInformation.stimOn)==0)
            stimData = neuralLFP(2:end,:);
        elseif(stimuli == 1) % all data before first stim
            stimData = neuralLFP(2:end,1:stimulationInformation.stimOn(stimuli)-1*30);
        elseif(stimuli == numel(stimulationInformation.stimOn)+1) % all data after last stim artifact
            stimData = neuralLFP(2:end,stimulationInformation.stimOn(stimuli-1)+5*30:end);
        elseif(stimulationInformation.stimOn(stimuli) - 5*30 > stimulationInformation.stimOn(stimuli-1)) % only count the data as non stim if 5 ms away from previous artifact
            stimData = neuralLFP(2:end,stimulationInformation.stimOn(stimuli-1)+5*30:stimulationInformation.stimOn(stimuli)-2*30);
        else
            stimData =[];
        end
        
        if(~isempty(stimData))
            try
                stimDataFiltered = acausalFilter(stimData')'; % filter the data 
                
                % compute threshold related information
                thresholdAll = thresholdAll + sum(stimDataFiltered.^2,2)/numPoints; % threshold based on SS data
                numPointsActual = size(stimDataFiltered,2) + numPointsActual;
                
            catch
                warning('thresholding error')
            end
        end
    end
    
    thresholdAll = sqrt(thresholdAll*numPoints/numPointsActual);
        
    %% build templates for each channel and waveform
    if(isfield('inputData','templateSubtract') && inputData.templateSubtract)
        templates = zeros(size(neuralLFP,1)-1,numel(unique(waveforms.waveSent)),size(neuralLFP,2));
        for stimuli = 2:numel(stimulationInformation.stimOn)
            stimData = neuralLFP(2:end,stimulationInformation.stimOn(stimuli-1):stimulationInformation.stimOn(stimuli)-10);
            stimDataFiltered = acausalFilter(stimData')'/sum(waveforms.waveSent(stimuli)==waveforms.waveSent);
            stimDataFiltered = reshape(stimDataFiltered,size(stimDataFiltered,1),1,size(stimDataFiltered,2));
            templates(:,waveforms.waveSent(stimuli),1:size(stimData,2)) = ...
                templates(:,waveforms.waveSent(stimuli),1:size(stimData,2)) + stimDataFiltered;
                
        end
    end
            
    %%

    spikeWaves = zeros(90000,lengthWave);
    spikeTimes = zeros(90000,1);
    spikeChan = zeros(90000,1);    

    spikeNum = 1;
    disp('extracting spikes')
    for ch = inputData.good_chan_list
        for stimIdx = 1:numel(stimulationInformation.stimOn)+1
            stimData = [];
            if(numel(stimulationInformation.stimOn) == 0)
                stimData = neuralLFP(ch+1,:);
            elseif(stimIdx == 1) % all data before first stim
                stimData = neuralLFP(ch+1,1:stimulationInformation.stimOn(stimIdx));
            elseif(stimIdx == numel(stimulationInformation.stimOn) + 1) % all data after last stimList
                stimData = neuralLFP(ch+1,stimulationInformation.stimOn(stimIdx-1):end);
            else % data before ith-1 stim up to ith stim
                stimData = neuralLFP(ch+1,stimulationInformation.stimOn(stimIdx-1):stimulationInformation.stimOn(stimIdx)-10);
            end

            % filter backwards on the channel
            stimData = acausalFilter(stimData');
            
            if(isfield('inputData','templateSubtract') && inputData.templateSubtract)
                stimData = stimData - squeeze(templates(ch,waveforms.waveSent(stimIdx),1:size(stimData,1)));
            end
            
            % compute threshold
            threshold = thresholdMult*thresholdAll(ch);
            if(abs(threshold) < 1) % if the threshold is super small, set it to a large value because the channel is bad
                threshold = sign(threshold)*100000;
            end

            %% deal with software amplification
            if(inputData.softwareAmplification && stimIdx > 1 && ch == inputData.dukeBoardChannel)
                stim_code = waveforms.waveSent(stimIdx-1);
                amp_idx = find(inputData.dukeBoard.amp_1_all == waveforms.parameters(stim_code).amp1);
                pw_idx = find(inputData.dukeBoard.pw1 == waveforms.parameters(stim_code).pWidth1); % in us
                polarity_idx = find(inputData.dukeBoard.polarity == waveforms.parameters(stim_code).polarity);
                
                wave_length = floor((waveforms.parameters(stim_code).pWidth1 + waveforms.parameters(stim_code).pWidth2 + ...
                    waveforms.parameters(stim_code).interphase)*(30000/10^6)); % in samples
                
                max_data_point_amplify = floor(max(inputData.dukeBoard.t_post_stim)*30000/10^3); % max_time_post_stim is in ms
                t_post_stim = (1:1:max_data_point_amplify)'/30;
                amplifier_gain = ones(size(t_post_stim));
                for t_idx = 1:numel(t_post_stim)
                    t_post_stim_idx = find(t_post_stim(t_idx) <= inputData.dukeBoard.t_post_stim,1,'first');
                    if(~isempty(t_post_stim_idx))
                        amplifier_gain(t_idx) = inputData.dukeBoard.gain_ratio(t_post_stim_idx,amp_idx,polarity_idx,pw_idx);
                    end
                end
                amplifier_gain(isnan(amplifier_gain)) = 1;
                
                stimData(wave_length+1:wave_length+max_data_point_amplify) = stimData(wave_length+1:wave_length+max_data_point_amplify)./amplifier_gain;
            end
            
            %% get threshold crossings
            thresholdCrossings = find(stimData>abs(threshold)); % acausal filtered data, positive threshold
            % check if too close to beginning or end
            if(isfield(inputData,'blankTime'))
                mask=thresholdCrossings>(preOffset+1+floor(inputData.blankTime*30)) & thresholdCrossings <numel(stimData)-(postOffset+1);
            else
                mask=thresholdCrossings>(preOffset+1) & thresholdCrossings <numel(stimData)-(postOffset+1);
            end
            
            thresholdCrossings = thresholdCrossings(mask);
            %% append data before and after stimData to get spikes near the edges
            numAppend = 100;
            stimData = [zeros(numAppend,size(stimData,2));stimData(:,:);zeros(numAppend,size(stimData,2))];
            thresholdCrossings = thresholdCrossings + numAppend;

            % remove potential artifacts based on max amplitude
            % of stim data
            % could be written to not be a for loop....
            crossingsMask = ones(numel(thresholdCrossings),1);
            for cross = 1:numel(thresholdCrossings)
                if(stimData(thresholdCrossings(cross),1) > maxAmplitude)
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
                    while crossCheck >= 1 && thresholdCrossings(crossCheck) >= thresholdCrossings(cross) - max(30,max(preOffset,postOffset))
                        crossingsMask(crossCheck) = 0;
                        crossCheck = crossCheck-1;
                    end
                end
            end  
            thresholdCrossings = thresholdCrossings(crossingsMask(:) == 1);

            % store thresholdCrossing data
            for cross = 1:numel(thresholdCrossings)
                if(stimIdx == 1)
                    spikeTimes(spikeNum) = neuralLFP(1,thresholdCrossings(cross)-numAppend); % this is in seconds
                else
                    spikeTimes(spikeNum) = neuralLFP(1,stimulationInformation.stimOn(stimIdx-1)+thresholdCrossings(cross)-1-numAppend); % this is in secondss
                end
                spikeChan(spikeNum) = ch;
                
                spikeWaves(spikeNum,:) = stimData((thresholdCrossings(cross)-preOffset):(thresholdCrossings(cross)+postOffset),1);
                spikeNum = spikeNum + 1;

                % preallocate space in arrays if close to max size
                if(spikeNum >= numel(spikeTimes))
                    spikeTimes = [spikeTimes;zeros(10000,1)];
                    spikeWaves = [spikeWaves; zeros(10000,lengthWave)];
                    spikeChan = [spikeChan; zeros(10000,1)];
                end
            end


        end % end channel for
    end % end stimOn for
        
    % store spike data
    disp('storing data')
    spikeTimes = spikeTimes(1:spikeNum-1,1);
    spikeWaves = spikeWaves(1:spikeNum-1,:);
    spikeChan = spikeChan(1:spikeNum-1,:);
    % sort the data by time
    [spikeTimes,sortOrder] = sort(spikeTimes);
    spikeWaves = spikeWaves(sortOrder,:);
    spikeChan = spikeChan(sortOrder,:);
    
    % store in nevData
    nevData.ts = spikeTimes;
    nevData.waveforms = spikeWaves(:,:);
    nevData.elec = spikeChan(:,:);
       
    % store raw data where a threshold crossing occurs (index of
    % spikeTimes in this file)
    rawData.ts = spikeTimes;
    rawData.waveforms = zeros(numel(spikeTimes),preOffset*2+postOffset*2+1);
    rawIdxList = zeros(numel(spikeTimes),1);
    stIdx = 1;
    for ttt = 1:size(neuralLFP,2)
        while(stIdx <= numel(spikeTimes) && neuralLFP(1,ttt) == spikeTimes(stIdx))
            rawIdxList(stIdx,1) = ttt;
            stIdx = stIdx + 1;
        end
    end
    for r = 1:numel(spikeTimes)
        rawIdx = rawIdxList(r);
        if(rawIdx ~= 0 && rawIdx-preOffset*2 > 0 && rawIdx+postOffset*2 <= size(neuralLFP,2))
            rawData.waveforms(r,:) = neuralLFP(spikeChan(r)+1,rawIdx-preOffset*2:rawIdx+postOffset*2)';
        end
    end
    rawData.elec = spikeChan;
        
    clear spikeTimes
    clear spikeChan
    clear spikeWaves
        
    % store artifact data
    artifactDataIndex = 1;
    artifactData.t = zeros(1000,1);
    % preset artifactData.artifact size based on num channels, subtract 1
    % for the time column
    artifactData.artifact = zeros(1000,numel(inputData.good_chan_list),artifactDataTime*30000/1000 + inputData.presample);
    
    stimulationInformation.stimOn = stimulationInformation.stimOn(stimOnMask==1);
    if(isfield(inputData,'artifactSkip'))
        artifactSkip = inputData.artifactSkip;
    else
        artifactSkip = 1;
        warning('storing artifact data for every stimulation which is very memory intensive');
    end
    
    for art = 1:artifactSkip:numel(stimulationInformation.stimOn)           
        if(stimulationInformation.stimOn(art)-inputData.presample > 1 && stimulationInformation.stimOn(art) + artifactDataTime*30000/1000 <= size(neuralLFP,2))
            artifactData.artifact(artifactDataIndex,:,:) = neuralLFP(inputData.good_chan_list+1,stimulationInformation.stimOn(art)-inputData.presample:stimulationInformation.stimOn(art)+floor(artifactDataTime*30000/1000)-1);
            artifactData.t(artifactDataIndex,1) = neuralLFP(1,stimulationInformation.stimOn(art));
            artifactDataIndex = artifactDataIndex + 1;
        end
        if(artifactDataIndex >= size(artifactData.artifact,1))
            artifactData.artifact(end+1:end+1000,:,:) = 0;
            artifactData.t(end+1:end+1000,1) = 0;
        end
    end
    
    artifactData.artifact = artifactData.artifact(1:artifactDataIndex-1,:,:);
    artifactData.t = artifactData.t(1:artifactDataIndex-1);
    
%     % store processed data for future use
%     save(strcat(inputData.folderpath,inputData.filename(1:end-4),'_nevData.mat'),'nevData','-v7.3');
%     save(strcat(inputData.folderpath,inputData.filename(1:end-4),'_rawData.mat'),'rawData','-v7.3');
%     save(strcat(inputData.folderpath,inputData.filename(1:end-4),'_artifactData.mat'),'artifactData','-v7.3');

    %% ADJUST STIM TIMES BASED ON RESETS IN DATA 
    % this is done for the spike times later, but the stim on data doesn't
    % change with it
    for resetIdx = 1:numel(outputData.DataPoints)
        stimulationInformation.stimOn(stimulationInformation.stimOn > outputData.DataPoints(resetIdx)) = ...
            stimulationInformation.stimOn(stimulationInformation.stimOn > outputData.DataPoints(resetIdx)) - outputData.DataPoints(resetIdx)-NSx_trim.MetaTags.Timestamp(resetIdx);
    end
    
    %% setup output data
    outputData.artifactData = artifactData;
    outputData.nevData = nevData;
    outputData.rawData = rawData;
    outputData.waveforms = waveforms;
    outputData.stimInfo = stimulationInformation;
    
end




% 
% %% fix stim times if joe sees an issue with the data file
%         if(inputData.issueExists)
%             % cable falling out, remove stim times that this corresponds
%             % too
%             maskArtifactKeep = ones(numel(stimulationInformation.stimOn),1);
%             for artIdx = 1:numel(stimulationInformation.stimOn)
%                 artData = cdsTempLFP(stimulationInformation.stimOn(artIdx):stimulationInformation.stimOn(artIdx)+90,inputData.dukeBoardChannel+1);
%                 if(mean(artData) < -8000 || max(artData-mean(artData)) < 500)
%                     maskArtifactKeep(artIdx) = 0;
%                 end
%             end
%             stimulationInformation.stimOn = stimulationInformation.stimOn(maskArtifactKeep==1);
%             % fix waveforms sent file
%             waveformFilename = strcat(fileList(i).name(1:underscoreIdx(end)),'waveformsSent',fname(underscoreIdx(end):end),'.mat');
%             load(waveformFilename)
%              
%             waveforms.chanSent = waveforms.chanSent(maskArtifactKeep == 1);
%             waveforms.waveSent = waveforms.waveSent(maskArtifactKeep == 1);
%             save(waveformFilename,'waveforms','-v7.3');
%         end