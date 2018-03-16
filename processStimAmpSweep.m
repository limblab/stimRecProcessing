function [ outputFigures,outputData ] = processStimAmpSweep(folderpath, inputData )
    %script to load stimulation files and generate perievent plots of 30khz
    %data. Formatted to work with runDataProcessing
    
    outputFigures=[];
    outputData=[];
    %get list of all files in the folder:
    if ~strcmp(folderpath(end),filesep)
        folderpath=[folderpath,filesep];
    end
    cd(folderpath);
    fileList=dir('*.nev');

    %get the mapfile info:
    %     mapFile=inputData.mapFile(8:end);
%     mapData=readtable(mapFile,'FileType','text','HeaderLines',13,'Delimiter','tab');
%     mapData.Properties.VariableNames{1}='col';%fix the column header
%     %fix the labels:
%     for i=1:size(mapData,1)
%         mapData.label{i}=[inputData.array1(6:end),mapData.label{i}];
%     end
    
    mapData=loadMapFile(inputData.mapFile(8:end));
        for i=1:size(mapData,1)
            mapData.label{i}=[inputData.array1(6:end),mapData.label{i}];
        end
    %establish trackers for information that is common to all files:
    eList={};
    posList=[];
    chList=[];
    if RDPIsAlreadyDone('artifactData',folderpath) && ~inputData.forceReload
        warning('processStimArtifact:foundExistingData','loading data from previous processing. This will have the PREVIOUS settings for time window, presample etc')
        artifactData =RDPLoadExisting('artifactData',folderpath);
        eList =RDPLoadExisting('eList',folderpath);
        posList =RDPLoadExisting('posList',folderpath);
        chList =RDPLoadExisting('chList',folderpath);
    else
        for i=1:numel(fileList)
            %% load file
            disp(['working on:'])
            disp(fileList(i).name)
            cds=commonDataStructure();
%             cds.file2cds([folderpath,fileList(i).name],inputData.ranBy,inputData.array1,inputData.monkey,inputData.lab,'ignoreJumps',inputData.task,inputData.mapFile,'recoverPreSync');
            cds.file2cds([folderpath,fileList(i).name],inputData.ranBy,inputData.array1,inputData.monkey,inputData.lab,'ignoreJumps',inputData.task,inputData.mapFile);
            %% find sync signal in analog data
            useSync=true;
            aIdx=[];
            syncIdx=[];
            if isempty(inputData.useSyncLabel)
            %look for sync under the label sync
                for j=1:numel(cds.analog)
                    syncIdx=find(strcmp(cds.analog{j}.Properties.VariableNames,'sync'));
                    if ~isempty(syncIdx)
                        aIdx=j;
                        syncName='sync';
                    end
                end
                %if it wasn't called sync, try for matt's 'StimTrig' label:
                if isempty(syncIdx)
                    for j=1:numel(cds.analog)
                        syncIdx=find(strcmp(cds.analog{j}.Properties.VariableNames,'StimTrig'));
                        if ~isempty(syncIdx)
                            aIdx=j;
                            syncName='StimTrig';
                        end
                    end
                end
                %if we didn't find a sync channel, just look for ainp16
                if isempty(syncIdx)
                    for j=1:numel(cds.analog)
                        syncIdx=find(strcmp(cds.analog{j}.Properties.VariableNames,'ainp16'));
                        if ~isempty(syncIdx)
                            useSync=false;
                            aIdx=j;
                            syncName='ainp16';
                        end
                    end
                end
                if isempty(aIdx)
                    error('processStimArtifact:cantFindSync','couldnt find a sync signal')
                end
            else
                useSync=inputData.useSyncLabel;
                if useSync
                    %find aIdx:
                    for j=1:numel(cds.analog)
                        syncIdx=find(strcmp(cds.analog{j}.Properties.VariableNames,'sync'));
                        if ~isempty(syncIdx)
                            aIdx=j;
                            syncName='sync';
                        end
                    end
                else
                    for j=1:numel(cds.analog)
                        syncIdx=find(strcmp(cds.analog{j}.Properties.VariableNames,'ainp16'));
                        if ~isempty(syncIdx)
                            useSync=false;
                            aIdx=j;
                            syncName='ainp16';
                        end
                    end
                end
            end
            %% use sync to get stim times:
            artifactData(i).stimOn=find(diff(cds.analog{aIdx}.(syncName)-mean(cds.analog{aIdx}.(syncName))>3)>.5);
%             if useSync
%                 stimOn=find(diff(cds.analog{aIdx}.sync-mean(cds.analog{aIdx}.sync)>3)>.5);
%             else
%                 stimOn=find(diff(cds.analog{aIdx}.ainp16-mean(cds.analog{aIdx}.ainp16)>3)>.5);
%             end
%             if numel(stimOn>100)
%                 %loop through till we find a likely channel:
%                 for i=2:size(cds.lfp,2)
%                 end
%             end
            stimOff=find(diff(cds.analog{aIdx}.(syncName)-mean(cds.analog{aIdx}.(syncName))<-3)>.5);
            artifactData(i).stimOff=nan(size(artifactData(i).stimOn));
            for j=1:numel(artifactData(i).stimOn)
                if j<numel(artifactData(i).stimOn)
                    next=artifactData(i).stimOn(j+1);
                else
                    next=numel(cds.analog{aIdx}.(syncName));
                end
                offIdx=stimOff(find((stimOff>artifactData(i).stimOn(j)& stimOff<next),1,'first'));
                if ~isempty(offIdx)
                    artifactData(i).stimOff(j)=offIdx;
                end
            end
            
            stimWindows=[artifactData(i).stimOn,artifactData(i).stimOn+inputData.windowSize-1];    
            %% find monitor data for cerestim monitoring ports:
            monitorName={'monitor','cerestim_module','Cerestim_module','Cerestim_Module'};
            mIdx=[];
            for j=1:numel(cds.analog)
                for k=1:numel(monitorName)
                    monitorIdx=find(strcmp(cds.analog{j}.Properties.VariableNames,monitorName{k}));
                    if ~isempty(monitorIdx)
                        mIdx=j;
                        monitorCol=monitorIdx;
                    end
                end
            end
            
            
            %% put spikes into structure:
            idxStart=strfind(fileList(i).name,'chan')+4;
            if numel(idxStart)>1
                %we have 2 instances of chan in the name, try again looking
                %for _chan
                idxStart=strfind(fileList(i).name,'_chan')+5;
            end
            idxEnd=strfind(fileList(i).name,'stim')-1;
            if numel(idxEnd)>1
                %we have 2 instances of end, look for stim_ instead and
                %hope that clears it up:
                idxEnd=strfind(fileList(i).name,'stim_')-1;
            end
            if isempty(idxStart)
                %hopefully this is a file that has the ch#stim format:
                for j=1:numel(idxEnd)
                    idxStart=strfind(fileList(i).name(1:idxEnd(j)),'ch');
                    if isempty(idxStart)
                        continue
                    end
                    idxStart=idxStart(end)+2;
                    if idxEnd(j)-idxStart<4
                        idxEnd=idxEnd(end);
                    end
                end
            end
            artifactData(i).stimChannel=str2num(fileList(i).name(idxStart:idxEnd));
            if isempty(artifactData(i).stimChannel)
                %we probably have the chan#_stim format:
                artifactData(i).stimChannel=str2num(fileList(i).name(idxStart:idxEnd-1));
            end

            if numel([artifactData.stimOff])<numel(inputData.minAmp:inputData.ampStep:inputData.maxAmp)
                error('failed to find the expected number of stim events. This will prevent proper amplitude assignment')
            end

%             labelList={cds.units.label};
%             arrayList={cds.units.array};
%             for j=1:numel(labelList)
%                 labelList{j}=[arrayList{j},labelList{j}];
%             end
            if isempty(cds.lfp)
                numChans=size(cds.analog{1,1},2);
            else
                numChans=size(cds.lfp,2);
            end
            artifactMat=[];
            monitorCh1Mat=[];
            for j=2:numChans
                for k=1:numel(artifactData(i).stimOn)
                    %disp(['j=',num2str(j),' k=',num2str(k)])
                    if isempty(cds.lfp)
                        artifactMat(j-1,k,1:inputData.windowSize+inputData.presample)=reshape(cds.analog{1,1}{stimWindows(k,1)-inputData.presample:stimWindows(k,2),j},[1,1,inputData.windowSize+inputData.presample]);
                        chanNum=str2num(cds.analog{1,1}.Properties.VariableNames{j}(5:end));
                        mapIdx=find(mapData.chan==chanNum);
                        electrodeList{j-1}=mapData.label{mapIdx};
                    else
                        artifactMat(j-1,k,1:inputData.windowSize+inputData.presample)=reshape(cds.lfp{stimWindows(k,1)-inputData.presample:stimWindows(k,2),j},[1,1,inputData.windowSize+inputData.presample]);
                        electrodeList{j-1}=cds.lfp.Properties.VariableNames{j};
                    end
                    %if we don't have a position for this electrode, find it:
                    if isempty(find(strcmp(eList,electrodeList{j-1}),1))
                        %find the electrode in the units structure and add it
                        %to eList and posList:
                        unitIdx=find(strcmp(mapData.label,electrodeList{j-1}));
                        if ~isempty(unitIdx)
                            eList(end+1)=electrodeList(j-1);
                            chList(end+1)=(uint8(mapData.bank{unitIdx})-65)*32+mapData.pin(unitIdx);
                            posList(end+1,:)=[mapData.row(unitIdx),mapData.col(unitIdx)];
                        end
                    end
                end
                %now get the monitor data showing the driving voltage the
                %stimulator used for the pulses
                if exist('mIdx','var') & ~isempty(mIdx)
                    for k=1:numel(artifactData(i).stimOn)
                        monitorCh1Mat(k,:)=cds.analog{mIdx}{(stimWindows(k,1)-inputData.presample):stimWindows(k,2)  , monitorCol}';
                    end
                end
            end
            artifactData(i).artifact=artifactMat;
            artifactData(i).electrodeNames=electrodeList;
            artifactData(i).monitor1=monitorCh1Mat;

            %clear cds
            clear cds
        end
        outputData.artifactData=artifactData;
        outputData.eList=eList;
        outputData.posList=posList;
        outputData.chList=chList;
    end
    
    %% calcualte stimulation current for each pulse
    baseSweep=[inputData.minAmp:inputData.ampStep:inputData.maxAmp]';
    numSweeps=numel(artifactData.stimOn)/numel(baseSweep);
    outputData.stimCurrent=repmat(baseSweep,[numSweeps,1]);
    %% plot artifacts grouped by stimulus current:
    chIdxMask=outputData.chList==outputData.stimChannel;
    for i=1:numel(baseSweep)
        %get indices of all stim events that match the current:
        StimIdxMask=outputData.stimCurrent==baseSweep(i);
        %plot all those data:
        outputFigures(end+1)=figure;
        plot(([inputData.presample:inputData.windowSize/30])/30,squeeze(outputData(chIdxMask,StimIdxMask,:))')
        title(['CH',num2str(outputData.stimChannel),' ',num2str(baseSweep(i)),'uA stim artifact'])
        set(outputFigures(end),'Name',['CH',num2str(outputData.stimChannel),'_',num2str(baseSweep(i)),'uA_stim_artifact'])
        xlabel('time after stim (ms)')
        ylabel('artifact amplitude (uV)')
    end
    
    
    %% calculate the rail time for each pulse: THIS IS BUGGY, DOES NOT ACCUNT FOR WHAT CURRENT WAS USED
%     for i=1:numel(artifactData)
%         [numElecs,numStimEvnts,pts]=size(artifactData(i).artifact);
%         artifactData(i).railTime=nan(numElecs,numStimEvnts);
%         artifactData(i).halfSettleTime=nan(numElecs,numStimEvnts);
%         for j=1:numElecs
%             for m=1:numStimEvnts
%                 railPoints=artifactData(i).artifact(j,m,:)<-8000|artifactData(i).artifact(j,m,:)>8000;
%                 midPoints=artifactData(i).artifact(j,m,:)<-4000|artifactData(i).artifact(j,m,:)>4000;
%                 railPoint=find(railPoints,1,'last');
%                 midPoint=find(midPoints,1,'last');
%                 if isempty(railPoint);
%                     railPoint=25;
%                 end
%                 if isempty(midPoint);
%                     midPoint=25;
%                 end
%                 if max(artifactData(i).artifact(j,m,railPoint:end))<8000
%                     %we Never hit the rail
%                     artifactData(i).railTime(j,m)=0;
%                     artifactData(i).halfSettleTime(j,m)=0;
%                 elseif min(artifactData(i).artifact(j,m,railPoint:end))>=8000
%                     %we stayed on the rail the whole time
%                     artifactData(i).railTime(j,m)=100;
%                     artifactData(i).halfSettleTime(j,m)=100;
%                 else
%                     artifactData(i).railTime(j,m)=railPoint/30;
%                     artifactData(i).halfSettleTime(j,m)=midPoint/30;
%                 end
%             end
%         end
%     end
%     %% get the stim voltage from the monitor pulse:
%     for i=1:numel(artifactData)
%         [numElecs,numStimEvnts,pts]=size(artifactData(i).artifact);    
%         artifactData(i).stimVoltage=nan(1,numStimEvnts);
%         artifactData(i).stimCurrent=reshape(repmat(1:100,[2,10]),[1,2000]);
%         warning('using only cathodal phase of stimulation to compute stim voltage')
%         for m=1:numStimEvnts
%             artifactData(i).stimVoltage(m)=max(artifactData(i).monitor1(m,:));
%         end
% 
%         % now find the voltage that is in the ADC range and create an
%         % extrapolated voltage:
%         inRange=artifactData(i).stimVoltage<4950 & artifactData(i).stimVoltage>-4950;
%         tmp=polyfit(artifactData(i).stimCurrent(inRange),artifactData(i).stimVoltage(inRange),1);
%         artifactData(i).voltageFactor=tmp(1);
%         artifactData(i).extrapVoltage=artifactData(i).stimCurrent*artifactData(i).voltageFactor;
%     end
    
    
    %% plots:
    if isfield(inputData,'doFigures') && ~inputData.doFigures
        return
    end
%     for i=1:numel(artifactData)
%         %railTime vs voltage
%         outputFigures(end+1)=figure;
%         voltageMat=reshape(artifactData(i).extrapVoltage(1:2:end),[100,10]);
%         railTimeMat=reshape(artifactData(i).railTime(artifactData(i).stimChannel,1:2:end),[100,10]);
%         plot(voltageMat,railTimeMat,'-*');
%         title(['rail time vs voltage for stim on CH',num2str(artifactData(i).stimChannel)])
%         set(outputFigures(end),'Name',['railTimeVSVoltage_stimCH',num2str(artifactData(i).stimChannel)])
%         xlabel('stim voltage (mV)')
%         ylabel('railTime (ms)')
%         %halfTime-railTime vs voltage
%         outputFigures(end+1)=figure;
%         settleTimeMat=reshape(artifactData(i).halfSettleTime(artifactData(i).stimChannel,1:2:end),[100,10]);
%         plot(voltageMat,settleTimeMat-railTimeMat,'-*');
%         title(['settleTime vs voltage for stim on CH',num2str(artifactData(i).stimChannel)]);
%         set(outputFigures(end),'Name',['settleTimeVSVoltage_stimCH',num2str(artifactData(i).stimChannel)]);
%         xlabel('stim voltage (mV)')
%         ylabel('time to settle to half rail from full rail (ms)');
% 
%     end
            


end

