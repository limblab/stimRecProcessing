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
    if(inputData.doAcausalFilter)
        disp('acausal high pass filtering')
        error('the following code has not been sanity checked for the correct transposition of data fed to the filter. do that check and remove this error')
        NSx.Data=acausalFilter(NSx.Data');
    else
        disp('high pass filtering')
        if median([NSx.ElectrodesInfo.HighFreqCorner])/1000<cutoff
            NSx.Data=filtfilt(b,a,NSx.Data')';
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
    %Data field
        %SerialDigitalIO
            nev.Data.SerialDigitalIO.InputType=[];
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
        %per nev file spec, header offset is number of bytes in header
        %basicHeaderBytes+extendedHeaderBytes
        %basic header bytes are 336
        %extended header bytes are:
        %32*(
        %       3*numChansWithSpikes    +  
        %       NumDigitalLabels        +  
        %       NumVideoSources         +
        %       NumTrackableObjectIDs   )
        %
        %there should always be 2 digital labels, 'serial' and 'digin' so
        %NumDigitalLabels will be 2
        %for files with no tracking or video, this condenses to:
        %336+32*(3*numChanWithSpikes+2)
        nev.MetaTags.HeaderOffset=336+32*(3*numel(nev.ElectrodesInfo)+numel(nev.IOLabels));
        nev.MetaTags.DataDuration=round(NSx.MetaTags.DataDurationSec*30000);
        nev.MetaTags.DataDurationSec=NSx.MetaTags.DataDurationSec;
        nev.MetaTags.PacketCount=1000;%this is a random value, I don't *think* it is saved, since openNEV computes it on the fly
        nev.MetaTags.TimeRes=sampleRate;
        nev.MetaTags.Application='File Dialog v6.03.01.00';
        nev.MetaTags.Filename=NSx.MetaTags.Filename;
        nev.MetaTags.FilePath=NSx.MetaTags.FilePath;
        nev.MetaTags.FileExt='.nev';
        nev.MetaTags.ChannelID=NSx.MetaTags.ChannelID;
    
    saveNEV(nev,[folderPath,'Output_Data',filesep, inputData.fileName(1:end-3),'nev'])
    outputData.nev=nev;
    outputData.ns5=NSx;
    figureList=[];