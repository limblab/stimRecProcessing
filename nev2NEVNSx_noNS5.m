function nev2NEVNSx_noNS5(fname,varargin)
    %this is a method function for the commonDataStructure class and should
    %be saved in the @commonDataStructure folder
    %
    %NEVNSx=nev2NEVNSx(fname)
    %loads data from the nev specified by the path in fname, and from
    %associated .nsx files with the same name. This code is derived from
    %cerebus2NEVNSx, but will NOT merge files, or look for keywords to pick
    %sorted files when the path to the unsorted file is given. fname must
    %be the FULL path and file name including extension
    %
    %this method is inteneded to be used internally during cds object
    %initiation, not called to generate NEVNSx objects in general.
    
    recoverPreSync=false;
    block='last';
    if ~isempty(varargin)
        for i=1:numel(varargin)
            if islogical(varargin{i})
                recoverPreSync=varargin{1};
            elseif ischar(varargin{i}) 
                switch varargin{i}
                    case 'first'
                        block='first';
                    case 'last'
                        block='last';
                    otherwise 
                        error('nev2NEVNSx:badBlockName','the block to use must be either first or last')
                end
            end
        end
        
        if recoverPreSync
            warning('nev2NEVNSx:dataRecoveryProbablyBad','the recover pre-sync data option is intended for files where the continuity of the data is not important (e.g. PESTH around an event). The data will be recovered by assuming a 100ms latency between data points. Analog data will be shifted and interpolated to fill missing points')
        end
    else
        recoverPreSync=false;
    end
    
    [folderPath,fileName,~]=fileparts(fname);
    
    %get the path for files matching the filename
    NEVpath = dir([folderPath filesep fileName '*.nev']);
    NSxList{1} = dir([folderPath filesep fileName '.ns1']);
    NSxList{2} = dir([folderPath filesep fileName '.ns2']);
    NSxList{3} = dir([folderPath filesep fileName '.ns3']);
    NSxList{4} = dir([folderPath filesep fileName '.ns4']);
    NSxList{5} = [];
    frequencies=[500 1000 2000 10000 30000];%vector of frequencies in the order that the NSx entries appear in NSxList
    
    %% populate cds.NEV
    if isempty(NEVpath)
        error('nev2NEVNSx:fileNotFound',['did not find a file with the path: ' fname])
    else
        if numel(NEVpath)>1
            %check to see if we have a sorted file with no digital:
            NEVpath=dir([folderPath filesep fileName '_nodigital*.nev']);
            digitalPath=dir([folderPath filesep fileName '_nospikes.mat']);
            sortedPath=dir([folderPath filesep fileName '-s.nev']);
            if ~isempty(NEVpath) && ~isempty(digitalPath)
                    if numel(NEVpath)>1
                        warning('nev2NEVNSx:multipleNodigitalFiles','found multiple files matching the *_nodigital.nev format. Attempting to identify the correct file')
                        nameLengths=cellfun(@length,{NEVpath.name});
                        NEVpath=NEVpath(nameLengths==max(nameLengths));
                        %now try to extract a number from the end of the
                        %path as we would see from the automatic
                        %save-scheme from plexon's offline sorter:
                        for i=1:numel(NEVpath)
                            fileNum=str2num(NEVpath(i).name(end-5:end-4));
                            if ~isempty(fileNum)
                                fileNumList(i)=fileNum;
                            else
                                fileNumList(i)=-10000;
                            end
                        end
                        NEVpath=NEVpath(find(fileNumList==max(fileNumList)));
                        disp(['continuing using file: ',NEVpath.name])
                    end
                    spikeNEV=openNEV('read', [folderPath filesep NEVpath.name],'nosave');
                    oldNEV=load([folderPath filesep digitalPath.name]);
                    oldNEVName=fieldnames(oldNEV);
                    oldNEV.(oldNEVName{1}).Data.Spikes=spikeNEV.Data.Spikes;
                    set(cds,'NEV',oldNEV.(oldNEVName{1}));
            elseif ~isempty(sortedPath)
                if numel(sortedPath)>1
                    error('nev2NEVNSx:multipleSorted',['found multiple sorted files in the target directory. Please remove the extraneous sorts, or rename them so that only 1 file has the format: FILENAME-s.nev'])
                end
                %check to see if we have a *.mat file with the same name as
                %our target:
                matPath=dir([folderPath filesep fileName '-s.mat']);
                if ~isempty(matPath)
                    disp(['located a mat-file with the same name as sorted file. Continuing using: ' folderPath filesep fileName '-s.mat'])
                    NEVpath=sortedPath;
                else
                    disp(['located a sorted file. Continuing using: ' folderPath filesep fileName '-s.nev'])
                    NEVpath=sortedPath;
                end
            else
                warning('nev2NEVNSx:multipleNEVFiles',['Found multiple files that start with the name given, but could not find files matching the pattern: ',fname,'_nodigital*.nev + ',fname,'_nospikes.mat'])
                disp(['continuing by loading the NEV that is an exact match for: ',fname,'.nev'])
                NEVpath = dir([folderPath filesep fileName '.nev']);
            end
%         else
%             set(cds,'NEV',openNEV('read', [folderPath filesep NEVpath.name],'nosave'));
        end
    end
    if ~exist('spikeNEV','var')
        %if we didn't load the NEV specially to merge digital data, load
        %the nev directly into the cds:
        set(cds,'NEV',openNEV('read', [folderPath filesep NEVpath.name],'nosave'));
    else
        
    end
    %identify resets in neural data:
    if ~isempty(cds.NEV)
        pData.numResets=numel(find(diff(double(cds.NEV.Data.Spikes.TimeStamp))<0));
        if pData.numResets>1
            warning('nev2NEVNSx:multipleResets','Multiple resync events found. This indicates a problem with the data file, please inspect manually.')
            disp('processing will continue assuming only the data after the last resync is valid')
            pData.resetTime=double(cds.NEV.Data.Spikes.TimeStamp(pData.numResets));
        elseif pData.numResets==1
            pData.resetTime=double(cds.NEV.Data.Spikes.TimeStamp(pData.numResets));
        end
        if pData.numResets>0
            syncIdxSpikes=find(diff(double(cds.NEV.Data.Spikes.TimeStamp))<0);
            syncIdxDigital=find(diff(double(cds.NEV.Data.SerialDigitalIO.TimeStamp))<0);
            if recoverPreSync
                for i=1:numel(syncIdxSpikes)
                    pData.stampShift(i)=double(cds.NEV.Data.Spikes.TimeStamp(syncIdxSpikes(i)))+round(.1/double(cds.NEV.MetaTags.SampleRes));
                    cds.NEV.Data.Spikes.TimeStamp(syncIdxSpikes(i)+1:end)=cds.NEV.Data.Spikes.TimeStamp(syncIdxSpikes(i)+1:end)+pData.stampShift(i);
                    %
                    if ~isempty(cds.NEV.Data.SerialDigitalIO.TimeStamp) && ~isempty(syncIdxDigital)%sometimes the sync will happen early enough
                        pData.timeShift(i)=double(cds.NEV.Data.SerialDigitalIO.TimeStampSec(syncIdxDigital(i)))+.1;
                        cds.NEV.Data.SerialDigitalIO.TimeStampSec(syncIdxDigital(i)+1:end)=cds.NEV.Data.SerialDigitalIO.TimeStampSec(syncIdxDigital(i)+1:end)+pData.timeShift(i);
                        cds.NEV.Data.Spikes.TimeStamp(syncIdxSpikes(i)+1:end)=cds.NEV.Data.Spikes.TimeStamp(syncIdxSpikes(i)+1:end)+pData.stampShift(i);
                    end
                end
            else
                %remove all timestamps before the sync events:
                spikeMask=false(numel(cds.NEV.Data.Spikes.TimeStamp),1);
                digitalMask=false(numel(cds.NEV.Data.SerialDigitalIO.TimeStampSec),1);
                
                for i=1:numel(syncIdxSpikes)
                    if strcmp(block,'last')%isolate data after last sync event
                        %spike timestamps:
                        spikeMask(1:syncIdxSpikes(i))=true;
                        %digital timestamps:
                        if ~isempty(cds.NEV.Data.SerialDigitalIO.TimeStamp) && ~isempty(syncIdxDigital)
                            digitalMask(1:syncIdxDigital(i))=true;
                        end
                    elseif strcmp(block,'first')%isolate data before first sync event
                        %spike timestamps:
                        spikeMask(syncIdxSpikes(i):end)=true;
                        %digital timestamps:
                        if ~isempty(cds.NEV.Data.SerialDigitalIO.TimeStamp) && ~isempty(syncIdxDigital)
                            digitalMask(syncIdxDigital(i):end)=true;
                        end
                    end
                end
                
                cds.NEV.Data.Spikes.TimeStamp(spikeMask)=[];
                cds.NEV.Data.Spikes.Unit(spikeMask)=[];
                cds.NEV.Data.Spikes.Electrode(spikeMask)=[];
                cds.NEV.Data.Spikes.Waveform(:,spikeMask)=[];
                cds.NEV.Data.SerialDigitalIO.TimeStamp(digitalMask)=[];
                cds.NEV.Data.SerialDigitalIO.TimeStampSec(digitalMask)=[];
                cds.NEV.Data.SerialDigitalIO.InsertionReason(digitalMask)=[];
                cds.NEV.Data.SerialDigitalIO.UnparsedData(digitalMask)=[];
                
            end
            cds.addProblem('detected reset events in the NEV where the cerebus clock reset to zero. ',pData)
        end
        
    end
    %% populate the cds.NSx fields
    for i=1:length(NSxList)
        fieldName=['NS',num2str(i)];
        if ~isempty(NSxList{i})
            %load the NSx into a temporary variable:
                %NSx=openNSxLimblab('read', [folderPath filesep NSxList{i}.name],'precision','short');
            NSx=openNSx('read', [folderPath filesep NSxList{i}.name],'precision','double','uv');
            %handle any dropped packets (cause data to show up as multiple
            %cells). sync with 2nd cerebus will also cause multiple cells
            %as output as the pre-sync and post-sync data go in different
            %cells
            if length(NSx.MetaTags.Timestamp)>1
                numResync=length(NSx.MetaTags.Timestamp)-1;
                %can't recall why I did the following, it seems extraneous but would skip very short data series:
%                 for idx = 1:length(NSx.MetaTags.Timestamp)-1
%                     if NSx.MetaTags.Timestamp(idx)+NSx.MetaTags.DataPoints(idx)*NSx.MetaTags.SamplingFreq>NSx.MetaTags.Timestamp(idx+1)
%                         numResync = numResync+1;
%                     end
%                 end
                if numResync>1
                    warning('nev2NEVNSx:multipleResets','Multiple resync events found. This indicates a problem with the data file, please inspect manually.')
                    disp(['continuing assuming only the data in the last cell of the NS',num2str(i),' is valid'])
                end
                if numResync<numel(NSx.MetaTags.Timestamp)-1
                    disp('This file may have packet loss. This file has less resync events than output cells.')
                end
                
                %add a note to the problems:
                pData.resetTimes=NSx.MetaTags.DataDurationSec(1);
                cds.addProblem(['detected reset events in the ',NSxList{i}.name,' where the cerebus clock reset to zero. Data is concatenated and there may be discontinuities'],pData)
                if isempty(pData.numResets) || numel(NSx.MetaTags.Timestamp)-1~=pData.numResets
                    warning('nev2NEVNSx:resetMismatch','the nev and the NSx have different numbers of time resets')
                    disp('this is PROBABLY due to sorting the nev and removing all the pre-sync data')
                    disp('cds loading will continue by eliminating pre-sync NSx data')
                    cds.addProblem(['reset mismatch between nev and ',NSxList{i}.name],NSx.MetaTags)
                end
                if recoverPreSync
                    DataDurationSec=NSx.MetaTags.DataDurationSec(1);
                    DataPointsSec=NSx.MetaTags.DataPointsSec(1);
                    DataPoints=NSx.MetaTags.DataPointsSec(1);
                    Data=NSx.Data{1};
                    for j=2:numResync+1
                        DataDurationSec=DataDurationSec+NSx.MetaTags.DataDurationSec(j)+.1;
                        DataPointsSec=DataPointsSec+NSx.MetaTags.DataPointsSec(j)+.1;
                        %get data for jth resync:
                        tmp=NSx.Data{j};
                        %construct timeseries for this resync:
                        offsetAdjust=NSx.MetaTags.Timestamp(j)-NSx.MetaTags.Timestamp(1);
                        tmpTime=offsetAdjust/30000+.1+[1:size(tmp,2)]/frequencies(i);
                        tgtTime=[1:1/frequencies(i):tmpTime(end)];
                        %interpolate data onto timeframe of first
                        %dataseries:
                        tmpData=interp1([0,tmpTime]',[Data(:,1),tmp]',tgtTime');
                        Data=[Data,tmpData'];
                    end
                    NSx.MetaTags.DataDurationSec=DataDurationSec;
                    NSx.MetaTags.DataPointsSec=DataPointsSec;
                    NSx.MetaTags.Timestamp=NSx.MetaTags.Timestamp(1);
                    NSx.Data=Data;
                    NSx.MetaTags.DataPoints=size(Data,2);
                else
                    if strcmp(block,'first')
                        NSx.MetaTags.DataDurationSec=NSx.MetaTags.DataDurationSec(1);
                        NSx.MetaTags.DataPoints=NSx.MetaTags.DataPoints(1);
                        NSx.MetaTags.DataPointsSec=NSx.MetaTags.DataPointsSec(1);
                        NSx.MetaTags.Timestamp=NSx.MetaTags.Timestamp(1);
                        NSx.Data=NSx.Data{1};
                    elseif strcmp(block,'last')
                        NSx.MetaTags.DataDurationSec=NSx.MetaTags.DataDurationSec(end);
                        NSx.MetaTags.DataPoints=NSx.MetaTags.DataPoints(end);
                        NSx.MetaTags.DataPointsSec=NSx.MetaTags.DataPointsSec(end);
                        NSx.MetaTags.Timestamp=NSx.MetaTags.Timestamp(end);
                        NSx.Data=NSx.Data{end};
                    end
                end
            end
                                    
            %insert into the cds
            set(cds,upper(fieldName),NSx)
        else
            %set the NSx field empty in case we are currently loading a
            % second NEV. This prevents re-loading data that was in one 
            %*.nev file but not the other when NEVNSx2cds is called
            set(cds,upper(fieldName),[])
        end
    end
    
    %%   now get info we will need to parse the NEVNSx data:
    NSxInfo.NSx_labels = {};
    NSxInfo.NSx_sampling = [];
    NSxInfo.NSx_idx = [];
    if ~isempty(cds.NS1)
        NSxInfo.NSx_labels = {NSxInfo.NSx_labels{:} cds.NS1.ElectrodesInfo.Label}';
        NSxInfo.NSx_sampling = [NSxInfo.NSx_sampling repmat(500,1,size(cds.NS1.ElectrodesInfo,2))];
        NSxInfo.NSx_idx = [NSxInfo.NSx_idx 1:size(cds.NS1.ElectrodesInfo,2)];
    end
    if ~isempty(cds.NS2)
        NSxInfo.NSx_labels = {NSxInfo.NSx_labels{:} cds.NS2.ElectrodesInfo.Label}';
        NSxInfo.NSx_sampling = [NSxInfo.NSx_sampling repmat(1000,1,size(cds.NS2.ElectrodesInfo,2))];
        NSxInfo.NSx_idx = [NSxInfo.NSx_idx 1:size(cds.NS2.ElectrodesInfo,2)];
    end
    if ~isempty(cds.NS3)
        NSxInfo.NSx_labels = {NSxInfo.NSx_labels{:} cds.NS3.ElectrodesInfo.Label};
        NSxInfo.NSx_sampling = [NSxInfo.NSx_sampling repmat(2000,1,size(cds.NS3.ElectrodesInfo,2))];
        NSxInfo.NSx_idx = [NSxInfo.NSx_idx 1:size(cds.NS3.ElectrodesInfo,2)];
    end
    if ~isempty(cds.NS4)
        NSxInfo.NSx_labels = {NSxInfo.NSx_labels{:} cds.NS4.ElectrodesInfo.Label}';
        NSxInfo.NSx_sampling = [NSxInfo.NSx_sampling repmat(10000,1,size(cds.NS4.ElectrodesInfo,2))];
        NSxInfo.NSx_idx = [NSxInfo.NSx_idx 1:size(cds.NS4.ElectrodesInfo,2)];
    end
    if ~isempty(cds.NS5)
        NSxInfo.NSx_labels = {NSxInfo.NSx_labels{:} cds.NS5.ElectrodesInfo.Label}';
        NSxInfo.NSx_sampling = [NSxInfo.NSx_sampling repmat(30000,1,size(cds.NS5.ElectrodesInfo,2))];
        NSxInfo.NSx_idx = [NSxInfo.NSx_idx 1:size(cds.NS5.ElectrodesInfo,2)];
    end
    %sanitize labels
    NSxInfo.NSx_labels = NSxInfo.NSx_labels(~cellfun('isempty',NSxInfo.NSx_labels));
    NSxInfo.NSx_labels = deblank(NSxInfo.NSx_labels);
    %apply aliases to labels:
    if ~isempty(cds.aliasList)
        for i=1:size(cds.aliasList,1)
            NSxInfo.NSx_labels(~cellfun('isempty',strfind(NSxInfo.NSx_labels,cds.aliasList{i,1})))=cds.aliasList(i,2);
        end
    end
    % check that we don't have a data stream using the reserved name
    % 'good'
    if ~isempty(find(strcmp('good',NSxInfo.NSx_labels),1));
        error('NEVNSx2cds:goodIsAReservedName','the cds and experiment code uses the label good as a flag for kinematic data, and treats this label specially when refiltering. This label is reserved to avoid unintended behaviro when refiltering other data sreams. Please use the alias function to re-name the good channel of input data')
    end
    
    set(cds,'NSxInfo',NSxInfo)
    
end