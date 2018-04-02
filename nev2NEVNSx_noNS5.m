function [outputData] = nev2NEVNSx_noNS5(fname,varargin)
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
    outputData = [];
    
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
                    outputData.oldNEV = oldNEV.(oldNEVName{1});
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
        outputData.NEV = openNEV('read', [folderPath filesep NEVpath.name],'nosave');
    else
        
    end
    %identify resets in neural data:
    if ~isempty(outputData.NEV)
        pData.numResets=numel(find(diff(double(outputData.NEV.Data.Spikes.TimeStamp))<0));
        if pData.numResets>1
            warning('nev2NEVNSx:multipleResets','Multiple resync events found. This indicates a problem with the data file, please inspect manually.')
            disp('processing will continue assuming only the data after the last resync is valid')
            pData.resetTime=double(outputData.NEV.Data.Spikes.TimeStamp(pData.numResets));
        elseif pData.numResets==1
            pData.resetTime=double(outputData.NEV.Data.Spikes.TimeStamp(pData.numResets));
        end
        if pData.numResets>0
            syncIdxSpikes=find(diff(double(outputData.NEV.Data.Spikes.TimeStamp))<0);
            syncIdxDigital=find(diff(double(outputData.NEV.Data.SerialDigitalIO.TimeStamp))<0);
            if recoverPreSync
                for i=1:numel(syncIdxSpikes)
                    pData.stampShift(i)=double(outputData.NEV.Data.Spikes.TimeStamp(syncIdxSpikes(i)))+round(.1/double(outputData.NEV.MetaTags.SampleRes));
                    outputData.NEV.Data.Spikes.TimeStamp(syncIdxSpikes(i)+1:end)=outputData.NEV.Data.Spikes.TimeStamp(syncIdxSpikes(i)+1:end)+pData.stampShift(i);
                    %
                    if ~isempty(outputData.NEV.Data.SerialDigitalIO.TimeStamp) && ~isempty(syncIdxDigital)%sometimes the sync will happen early enough
                        pData.timeShift(i)=double(outputData.NEV.Data.SerialDigitalIO.TimeStampSec(syncIdxDigital(i)))+.1;
                        outputData.NEV.Data.SerialDigitalIO.TimeStampSec(syncIdxDigital(i)+1:end)=outputData.NEV.Data.SerialDigitalIO.TimeStampSec(syncIdxDigital(i)+1:end)+pData.timeShift(i);
                        outputData.NEV.Data.Spikes.TimeStamp(syncIdxSpikes(i)+1:end)=outputData.NEV.Data.Spikes.TimeStamp(syncIdxSpikes(i)+1:end)+pData.stampShift(i);
                    end
                end
            else
                %remove all timestamps before the sync events:
                spikeMask=false(numel(outputData.NEV.Data.Spikes.TimeStamp),1);
                digitalMask=false(numel(outputData.NEV.Data.SerialDigitalIO.TimeStampSec),1);
                
                for i=1:numel(syncIdxSpikes)
                    if strcmp(block,'last')%isolate data after last sync event
                        %spike timestamps:
                        spikeMask(1:syncIdxSpikes(i))=true;
                        %digital timestamps:
                        if ~isempty(outputData.NEV.Data.SerialDigitalIO.TimeStamp) && ~isempty(syncIdxDigital)
                            digitalMask(1:syncIdxDigital(i))=true;
                        end
                    elseif strcmp(block,'first')%isolate data before first sync event
                        %spike timestamps:
                        spikeMask(syncIdxSpikes(i):end)=true;
                        %digital timestamps:
                        if ~isempty(outputData.NEV.Data.SerialDigitalIO.TimeStamp) && ~isempty(syncIdxDigital)
                            digitalMask(syncIdxDigital(i):end)=true;
                        end
                    end
                end
                
                outputData.NEV.Data.Spikes.TimeStamp(spikeMask)=[];
                outputData.NEV.Data.Spikes.Unit(spikeMask)=[];
                outputData.NEV.Data.Spikes.Electrode(spikeMask)=[];
                outputData.NEV.Data.Spikes.Waveform(:,spikeMask)=[];
                outputData.NEV.Data.SerialDigitalIO.TimeStamp(digitalMask)=[];
                outputData.NEV.Data.SerialDigitalIO.TimeStampSec(digitalMask)=[];
                outputData.NEV.Data.SerialDigitalIO.InsertionReason(digitalMask)=[];
                outputData.NEV.Data.SerialDigitalIO.UnparsedData(digitalMask)=[];
                
            end
        end
        
    end
    
  
    
    
end