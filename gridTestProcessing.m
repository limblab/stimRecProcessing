%% set up
targetDir='/media/tucker/My Passport/local processing/stimTesting/20180227_stimswitch/';


%% parse files into folders:
clear fileParams
    tmp=dir(targetDir);
    fileNames={tmp.name};
    
    for i=1:numel(fileNames)
        
        A1Loc=strfind(fileNames{i},'_A1-');
        if isempty(A1Loc)
            continue
        end
        A2Loc=strfind(fileNames{i},'_A2-');
        PW1Loc=strfind(fileNames{i},'_PW1-');
        PW2Loc=strfind(fileNames{i},'_PW2-');
        interpulseLoc(1)=strfind(fileNames{i},'_interpulse');
        tmp=strfind(fileNames{i}(interpulseLoc(1)+11:end),'_');
        interpulseLoc(2)=interpulseLoc(1)+11+tmp(1)-1;
        chanLoc(1)=strfind(fileNames{i},'_chan');
        tmp=strfind(fileNames{i}(chanLoc(1):end),'stim_');
        chanLoc(2)=chanLoc(1)+tmp(1)-1;
        %get A1
        fileParams(i).A1=str2num(fileNames{i}(A1Loc+4:A2Loc-1));
        %get A2
        fileParams(i).A2=str2num(fileNames{i}(A2Loc+4:PW1Loc-1));
        %get PW1
        fileParams(i).PW1=str2num(fileNames{i}(PW1Loc+5:PW2Loc-1));
        %get PW2
        fileParams(i).PW2=str2num(fileNames{i}(PW2Loc+5:interpulseLoc(1)-1));
        %get interpulse
        fileParams(i).interpulse=str2num(fileNames{i}(interpulseLoc(1)+11:interpulseLoc(2)-1));
        %get channel
        fileParams(i).chan=str2num(fileNames{i}(chanLoc(1)+5:chanLoc(2)-1));
        fileParams(i).name=fileNames{i};
    end
    %prune irrelevant entries:
    for i=numel(fileParams):-1:1
        if isempty(fileParams(i).A1)
            fileParams(:,i)=[];
        end
    end
    %create proper folders
    
    interpulseList=unique([fileParams.interpulse]);
    chanList=unique([fileParams.chan]);
    
    %move into appropriate folder
    for i=1:numel(fileParams)
        IP=['IP',num2str(fileParams(i).interpulse)];
        if isempty(dir([targetDir IP]))
            %make the interpulse top level directory
            mkdir(targetDir, IP)
        end
        chan=['CH',num2str(fileParams(i).chan)];
        if isempty(dir([targetDir IP filesep chan]))
            %make the chan second level directory
            mkdir([targetDir IP],chan);
        end
        %compose folder string:
        folderString=['A1-',num2str(fileParams(i).A1),...
                        '_A2-',num2str(fileParams(i).A2),...
                        '_PW1-',num2str(fileParams(i).PW1),...
                        '_PW2-',num2str(fileParams(i).PW2)];
        if isempty(dir([targetDir IP filesep chan filesep folderString]))
            %make the directory for this test
            mkdir([targetDir IP filesep chan],folderString)
        end
        %copy this test into the directory:
        movefile([targetDir fileParams(i).name], [targetDir IP filesep chan filesep folderString filesep fileParams(i).name])
    end
%% process tests of interest:
functionName='processStimArtifact';

inputData.task='tasknone';
inputData.ranBy='ranByTucker'; 
inputData.array1='arrayResistor'; 
inputData.monkey='monkeyResistor';
%han
% inputData.mapFile='mapFile/media/tucker/My Passport/local processing/stimTesting/20161112/SN 6251-001459.cmp';
%chips
% inputData.mapFile='mapFile/media/tucker/My Passport/local processing/stimTesting/SN 6251-001455.cmp';
%chewie
%inputData.mapFile='mapFile/media/tucker/My Passport/local processing/stimTesting/20161205_chewie_PMDStim_PMD-recording/Chewie Left PMd SN 6251-001469.cmp';
%saline
inputData.mapFile='mapFile/media/tucker/My Passport/local processing/stimTesting/20161020_saline/1025-0370.cmp';
%saline2
%inputData.mapFile='mapFile/media/tucker/My Passport/local processing/stimTesting/20161220_saline/SN 6251-001695.cmp';
inputData.badChList=[];
inputData.windowSize=30*30;%in points. multiply ms by 30 to get points
inputData.presample=5;%in points
inputData.plotRange=8.3;%in mV
inputData.interpulse=.000053;%in s
inputData.lab=6;
inputData.useSyncLabel=[];
inputData.doFilter=false;
inputData.syncLength=0;%.00000200;%in s
inputData.forceReload=true;    

%get list of interpulse
IPFolders=dir(targetDir);
IPList=unique({IPFolders.name});
%loop through interpulse
for i=3:numel(IPList)
    inputData.interpulse=str2num(IPList{i}(3:end))/1000000;
    %get list of stimulated channels:
    chanFolders=dir([targetDir IPList{i}]);
    chanList=unique({chanFolders.name});
    %loop through channels
    for j=3:min([numel(chanList),6])
        %get list of tests
        testFolders=dir([targetDir IPList{i} filesep chanList{j}]);
        testNames={testFolders.name};
        %loop through test folders
        for k=3:numel(testFolders)
            %parse test for parameters
            PW1Loc=strfind(testNames{k},'_PW1-');
            PW2Loc=strfind(testNames{k},'_PW2-');
            inputData.pWidth1=str2num(testNames{k}(PW1Loc+5:PW2Loc-1))/1000000;
            inputData.pWidth2=str2num(testNames{k}(PW2Loc+5:end))/1000000;
            %initiate runDataProcessing 
            dataStruct2 = runDataProcessing(functionName,[targetDir IPList{i} filesep chanList{j} filesep testNames{k}],inputData);
            close all
            clear dataStruct2
        end
    end
end


%%  make plots

t=[1:size(artifactData.artifact,3)]/30;
figure;subplot(6,1,1); for i=1:2:10; plot(t,squeeze(artifactData.artifact(1,i,:)),'r');hold on;end
subplot(6,1,1); for i=2:2:10; plot(t,squeeze(artifactData.artifact(1,i,:)),'b');hold on;end
subplot(6,1,2); for i=1:2:10; plot(t,squeeze(artifactData.artifact(2,i,:)),'r');hold on;end
subplot(6,1,2); for i=2:2:10; plot(t,squeeze(artifactData.artifact(2,i,:)),'b');hold on;end
subplot(6,1,3); for i=1:2:10; plot(t,squeeze(artifactData.artifact(3,i,:)),'r');hold on;end
subplot(6,1,3); for i=2:2:10; plot(t,squeeze(artifactData.artifact(3,i,:)),'b');hold on;end
subplot(6,1,4); for i=1:2:10; plot(t,squeeze(artifactData.artifact(4,i,:)),'r');hold on;end
subplot(6,1,4); for i=2:2:10; plot(t,squeeze(artifactData.artifact(4,i,:)),'b');hold on;end
subplot(6,1,5); for i=1:2:10; plot(t,squeeze(artifactData.artifact(5,i,:)),'r');hold on;end
subplot(6,1,5); for i=2:2:10; plot(t,squeeze(artifactData.artifact(5,i,:)),'b');hold on;end
subplot(6,1,6); for i=1:2:10; plot(t,squeeze(artifactData.artifact(6,i,:)),'r');hold on;end
subplot(6,1,6); for i=2:2:10; plot(t,squeeze(artifactData.artifact(6,i,:)),'b');hold on;end


% t=[0:30]/30;
% figure;subplot(6,1,1); for i=1:2:10; plot(t,squeeze(artifactData.artifact(1,i,5:35)),'r');hold on;end
% subplot(6,1,1); for i=2:2:10; plot(t,squeeze(artifactData.artifact(1,i,5:35)),'b');hold on;end
% subplot(6,1,2); for i=1:2:10; plot(t,squeeze(artifactData.artifact(2,i,5:35)),'r');hold on;end
% subplot(6,1,2); for i=2:2:10; plot(t,squeeze(artifactData.artifact(2,i,5:35)),'b');hold on;end
% subplot(6,1,3); for i=1:2:10; plot(t,squeeze(artifactData.artifact(3,i,5:35)),'r');hold on;end
% subplot(6,1,3); for i=2:2:10; plot(t,squeeze(artifactData.artifact(3,i,5:35)),'b');hold on;end
% subplot(6,1,4); for i=1:2:10; plot(t,squeeze(artifactData.artifact(4,i,5:35)),'r');hold on;end
% subplot(6,1,4); for i=2:2:10; plot(t,squeeze(artifactData.artifact(4,i,5:35)),'b');hold on;end
% subplot(6,1,5); for i=1:2:10; plot(t,squeeze(artifactData.artifact(5,i,5:35)),'r');hold on;end
% subplot(6,1,5); for i=2:2:10; plot(t,squeeze(artifactData.artifact(5,i,5:35)),'b');hold on;end
% subplot(6,1,6); for i=1:2:10; plot(t,squeeze(artifactData.artifact(6,i,5:35)),'r');hold on;end
% subplot(6,1,6); for i=2:2:10; plot(t,squeeze(artifactData.artifact(6,i,5:35)),'b');hold on;end
