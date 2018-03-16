function [figureList,outputData]=getSpikesFromFullBW(folderPath,inputData)
%% setup:
    %standard deviations from mean for theshold
    stdErrThresh=inputData.stdErrThresh;

    %clear high channels flag
    clearHighChans=inputData.clearHighChans;
    %snippitWindow
    preSample=inputData.preSample;
    postSample=inputData.postSample;
    snippit=[-preSample:postSample];
                
    sampleRate=30000;
    %file to work on
%     filename=
    %filter
    cutoff=inputData.HPFreq;%hz
    nPoles=inputData.poles;%1/2 the effective poles since filtfilt doubles effective number of poles
    [b,a]=butter(nPoles,2*cutoff/sampleRate,'high');
    
    
%% load file:
disp('opening file')
    NSx=openNSx('read', [folderPath,inputData.fileName]);
%% clear unwanted high channels:
if clearHighChans
    disp('removing unwanted channels')
    if max([NSx.ElectrodesInfo.ElectrodeID])>96
        chanMask=[NSx.ElectrodesInfo.ElectrodeID]<97;
        %data
            NSx.Data=NSx.Data(chanMask,:);
        %meta
            NSx.MetaTags.ChannelID=NSx.MetaTags.ChannelID(chanMask);
            NSx.MetaTags.ChannelCount=sum(chanMask);
        %electrode info
            NSx.ElectrodesInfo=NSx.ElectrodesInfo(chanMask);
    end
end
%% HP filter
    disp('high pass filtering')
    if median([NSx.ElectrodesInfo.HighFreqCorner])/1000<cutoff
        NSx.Data=filtfilt(b,a,NSx.Data')';
    end

%% select good channels:
if inputData.manualCheck
    chanMask=true(size(NSx.Data,1),1);    
    for i=1:size(NSx.Data,1)
        tmpData=reshape(NSx.Data(i,1:10000)/4,[100,100]);
        figure;plot(tmpData )
        inputChar='?';
        while ~strcmpi(inputChar,'g') && ~strcmpi(inputChar,'b')
            inputChar=input('is the channel good?','s');
            if numel(inputChar)>1
                warning('only single characters g or b accepted')
            end
        end
        if strcmpi(inputChar,'b')
            chanMask(i)=0;
        end
    end
end

%% rereference good channels to mean of good channels
if inputData.meanFilter
    disp('re-referencing to the common mode mean')
    NSx.Data=NSx.Data(chanMask,:)-repmat(mean(NSx.Data(chanMask,:)),[sum(chanMask),1]);
    %also clean up bad channels from meta and electrodsInfo
    NSx.MetaTags.ChannelID=NSx.MetaTags.ChannelID(chanMask);
    NSx.MetaTags.ChannelCount=sum(chanMask);
    NSx.ElectrodesInfo=NSx.ElectrodesInfo(chanMask);
end
%% PCA filter
if inputData.PCAFilter
    nPoints=size(NSx.Data,2);
    disp('estimating PCs of the data')
    [coeffData]=pca(single(NSx.Data(:,randsample(nPoints,floor(nPoints/10)))'));
    
    disp('shuffling the data')
    %create mask to shuffle data to generate baseline PC weights:
    mask=zeros(size(NSx.Data,1),floor(nPoints/10));
    for i=1:size(mask,1);
        mask(i,:)=datasample(1:nPoints,floor(nPoints/10),'Replace',false);%datasample appears to be *slightly* faster than randsample for large data
    end
    %get PCs of shuffled data
    disp('calculating PCs of shuffled data')
    coeffBase=pca(single(NSx.Data(mask)'));
    %compare PC values of real data to shuffled distribution 
    disp('checking whether weights in PCs of real data fall outside the normal range for shuffled data')
    baseDist=fitdist(coeffBase(:),'kernel');
    probs=reshape(cdf(baseDist,coeffData(:)),size(coeffData));
    mask=sum(probs>.95)>5;
    numPCs=sum(mask);
    if numPCs>0
        disp(['found ',num2str(numPCs),' PCs that have weights across multiple channels'])
        disp('removing these PCs')
        tmp=double(NSx.Data')*inv(coeffData);
        tmp(:,mask)=0;
        NSx.Data=int16(round(NSx.Data-int16(tmp*coeffData)'));
    end
end
%% threshold
disp('find threshold crossings')
    thresholdCrossings=nan(size(NSx.Data,1), ceil(size(NSx.Data,2)/postSample));
    for i=1:size(NSx.Data,1)
        thresholdVal(i)=-1*stdErrThresh*sqrt(var(single(NSx.Data(i,:))));
        if NSx.Data(i,1)<thresholdVal(i)
            descending=false;
        else
            descending=true;
        end
        n=1;
        for j=preSample+1:size(NSx.Data,2)-postSample-1
            if NSx.Data(i,j)<thresholdVal(i) && descending
                %log point in timeseries
                thresholdCrossings(i,n)=j;
                n=n+1;
                %skip ahead by sample window
                j=j+postSample+1;
                %set descending flag to false so we don't re-trigger on the
                %   same wave
                descending =false;
            elseif NSx.Data(i,j)>thresholdVal(i) && ~descending
                %reset descending flag
                descending=true;
            end
        end
    end

%% compose nev
disp('starting to compose nev structure')
    %MetaTags field
        nev.MetaTags.Subject=[];
        nev.MetaTags.Experimenter=[];
        nev.MetaTags.DateTime=NSx.MetaTags.DateTime;
        nev.MetaTags.SampleRes=uint32(NSx.MetaTags.SamplingFreq);
        nev.MetaTags.Comment=NSx.MetaTags.Comment;
        nev.MetaTags.FileTypeID='NEURALEV';
        nev.MetaTags.Flags='0000000000000001';
        nev.MetaTags.openNEVver=[];
        nev.MetaTags.DateTimeRaw=NSx.MetaTags.DateTimeRaw;
        nev.MetaTags.FileSpec='2.3';
        nev.MetaTags.HeaderOffset=14224;
        nev.MetaTags.DataDuration=round(NSx.MetaTags.DataDurationSec*30000);
        nev.MetaTags.DataDurationSec=NSx.MetaTags.DataDurationSec;
        nev.MetaTags.PacketCount=1000;%this is a random value, I don't *think* it is saved, since openNEV computes it on the fly
        nev.MetaTags.TimeRes=sampleRate;
        nev.MetaTags.Application='File Dialog v6.03.01.00';
        nev.MetaTags.Filename=NSx.MetaTags.Filename;
        nev.MetaTags.FilePath=NSx.MetaTags.FilePath;
        nev.MetaTags.FileExt='.nev';
        nev.MetaTags.ChannelID=NSx.MetaTags.ChannelID;
    %Data field
        %SerialDigitalIO
            nev.Data.SeriNSx.DataalDigitalIO.InputType=[];
            nev.Data.SerialDigitalIO.TimeStamp=[];
            nev.Data.SerialDigitalIO.TimeStampSec=[];
            nev.Data.SerialDigitalIO.Type=[];
            nev.Data.SerialDigitalIO.Value=[];
            nev.Data.SerialDigitalIO.InsertionReason=[];
            nev.Data.SerialDigitalIO.UnparsedData=[];
        %Spikes
            for i= 1:size(NSx.Data,1)
                points{i}=thresholdCrossings(i,~isnan(thresholdCrossings(i,:)));
                sampleMat=repmat(snippit,[numel(points{i}),1])+repmat(points{i}',[1,numel(snippit)]);
                currentWaves{i}=reshape(NSx.Data(i,sampleMat(:)),size(sampleMat));
                electrodeList{i}=NSx.ElectrodesInfo(i).ElectrodeID*uint16(ones(size(points{i})));
            end
            nev.MetaTags.PacketBytes=2*(size(currentWaves{1},2)+4);
            %TimeStamp
                nev.Data.Spikes.TimeStamp=uint32(cell2mat(points));
            %Electrode
                nev.Data.Spikes.Electrode=uint16(cell2mat(electrodeList));
            %Unit
                nev.Data.Spikes.Unit=uint8(zeros(size(nev.Data.Spikes.Electrode)));
            %Waveform
                nev.Data.Spikes.Waveform=int16(cell2mat(currentWaves'))';
            %WaveformUnit
                nev.Data.Spikes.WaveformUnit='raw';
            
        %Comments
            nev.Data.Comments.TimeStampStarted=[];
            nev.Data.Comments.TimeStampStartedSec=[];
            nev.Data.Comments.TimeStamp=[];
            nev.Data.Comments.TimeStampSec=[];
            nev.Data.Comments.CharSet=[];
            nev.Data.Comments.Text=[];
        %VideoSync
            nev.Data.VideoSync.TimeStamp=[];
            nev.Data.VideoSync.FileNumber=[];
            nev.Data.VideoSync.FrameNumber=[];
            nev.Data.VideoSync.ElapsedTime=[];
            nev.Data.VideoSync.SourceID=[];
        %Tracking
            nev.Data.Tracking=[];
        %TrackingEvents
            nev.Data.TrackingEvents.TimeStamp=[];
            nev.Data.TrackingEvents.TimeStampSec=[];
            nev.Data.TrackingEvents.ROIName=[];
            nev.Data.TrackingEvents.ROINum=[];
            nev.Data.TrackingEvents.Event=[];
            nev.Data.TrackingEvents.Frame=[];
        %PatientTrigger
            nev.Data.PatientTrigger.TimeStamp=[];
            nev.Data.PatientTrigger.TriggerType=[];
        %Reconfig
            nev.Data.Reconfig.TimeStamp=[];
            nev.Data.Reconfig.ChangeType=[];
            nev.Data.Reconfig.ChangeName=[];
            nev.Data.Reconfig.ConfigChanged=[];
    %ElectrodesInfo field
        for i=1:size(NSx.Data,1)
            nev.ElectrodesInfo(i).ElectrodeID=NSx.ElectrodesInfo(i).ElectrodeID;
            nev.ElectrodesInfo(i).ConnectorBank=NSx.ElectrodesInfo(i).ConnectorBank;
            nev.ElectrodesInfo(i).ConnectorPin=NSx.ElectrodesInfo(i).ConnectorPin;
            nev.ElectrodesInfo(i).DigitalFactor=250;
            nev.ElectrodesInfo(i).EnergyThreshold=0;
            nev.ElectrodesInfo(i).HighThreshold=0;
            nev.ElectrodesInfo(i).LowThreshold=thresholdVal(i);
            nev.ElectrodesInfo(i).Units=0;
            switch class(nev.Data.Spikes.Waveform)
                case 'int8'
                    nev.ElectrodesInfo(i).WaveformBytes=1;
                case 'int16'
                    nev.ElectrodesInfo(i).WaveformBytes=2;
                case 'int32'
                    nev.ElectrodesInfo(i).WaveformBytes=4;
                otherwise
                    error('data must be cast as an integer, typically int16')
            end
            nev.ElectrodesInfo(i).ElectrodeLabel=NSx.ElectrodesInfo(i).Label;
            nev.ElectrodesInfo(i).HighFreqCorner=NSx.ElectrodesInfo(i).HighFreqCorner;
            nev.ElectrodesInfo(i).HighFreqOrder=NSx.ElectrodesInfo(i).HighFreqOrder;
            nev.ElectrodesInfo(i).HighFilterType=NSx.ElectrodesInfo(i).HighFilterType;
            nev.ElectrodesInfo(i).LowFreqCorner=NSx.ElectrodesInfo(i).LowFreqCorner;
            nev.ElectrodesInfo(i).LowFreqOrder=NSx.ElectrodesInfo(i).LowFreqOrder;
            nev.ElectrodesInfo(i).LowFilterType=NSx.ElectrodesInfo(i).LowFilterType;
        end
    %IOLabels field
    nev.IOLabels={['serial' 0 0 0 0 0 0 0 0 0 0],['digin' 0 0 0 0 0 0 0 0 0 0 0] };
    
    saveNEV(nev,[folderPath,'Output_Data',filesep, inputData.fileName(1:end-3),'nev'])
    outputData.nev=nev;
    outputData.ns5=NSx;
    figureList=[];