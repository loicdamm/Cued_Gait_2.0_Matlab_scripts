%% COMPIL_CUED_GAIT_DATA 2023-2024
clear all

%% Paths
addpath CircStat2012a/

%% TO DOs
% universal import .mat .csv 30/01/24
% sampling rate 31/01/24
% output variables: 2 series / 0 01/02/24
% output variables: suppress 0 01/02/24
% Conversion variables in ms 01/02/24
% Ref cadence
% relative phase

%% Participants' directories

DirFilesList = dir('../DATA/');
DirFlags = [DirFilesList.isdir];
% Extract only those that are directories.
subFolders = DirFilesList(DirFlags); % A structure with extra info.
% Get only the folder names into a cell array.
%subFolderNames = {subFolders(3:end-1).name}; % Start at 3 to skip . and .. and exclude the last one "discarded"
subFolders=subFolders(contains({subFolders.name},'0'));
subFolderNames = {subFolders(:).name};

%% Import participants' Kq

part_Kq = read_part_Kq("../DATA/données participants.xlsx");
condiList={'Tapping', 'Gait'};
% new variable condiListFileName since tapping acquired with coda and Gait
% data acquired with delsys system
condiListFileName={'coda', 'delsys'};

%% Definition of constants
nbOfSteps=9;

%% PREALLOCATION
cuedGait(28,2,5).CT=nan(1,nbOfSteps);
cuedGait(28,2,5).ST=nan(1,nbOfSteps);
cuedGait(28,2,5).relPh=nan(1,nbOfSteps);

%% Participant loop

for iPart = 1 : length(subFolderNames)
    fprintf('\n Sub folder #%d = %s\n', iPart, subFolderNames{iPart});
    % Participant' files
    %filesList = dir(fullfile('../DATA/',subFolderNames{iPart}, '/*.mat'));
    filesList = dir(fullfile('../DATA/',subFolderNames{iPart}));
    nPart=find(contains(part_Kq{:,3},subFolderNames{iPart}));
    bpmTap=part_Kq{nPart,7};
    cad=part_Kq{nPart,10};
    refIoiWalk=60/cad*1000;
    refIoiTap=60/bpmTap*1000;
    %refIoiTap=part_Kq{nPart,9};


    %% Condition Loop 1=tapping 2 = gait

    for noCondi=1:2

        %% Files Loop Gait 3 files Tapping 4 files
        %fileInd=strfind({filesList.name}, condiList{noCondi});
        Desc   = {filesList.name};
        Desc(~cellfun('isclass', Desc, 'char')) = {''};
        matchC = reshape(strfind(Desc, condiListFileName{noCondi}), size(filesList));
        match  = ~cellfun('isempty', matchC);
        fileInd=find(match==1);
        if isempty(fileInd)==1
            fprintf('Warning!!! One condition missing: %s\n', condiList{noCondi})
        end

        %% Tapping
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if noCondi==1%tapping
            peakDAna=125;%round(60/bpmTap*250);

            for iFile=1:length(fileInd)
                %% Retrouver les dephasages du signal

                file2read=fullfile('../DATA/',subFolderNames{iPart},filesList(fileInd(iFile)).name);
                fprintf('Reading %s\n', file2read)

                load(file2read);
                %                 for iiii=1:30
                %
                %                     figure(iPart)
                %                     plot(Analog.Analog01.value)
                %                     hold on
                %                     plot(Analog.Analog02.value)
                %                     plot(Analog.Analog03.value)
                %                     plot(Analog.Analog04.value)
                %                     plot(Analog.Analog05.value)
                %                     plot(Analog.Analog06.value)
                %                     plot(Analog.Analog07.value)
                %                     plot(Analog.Analog08.value)
                %                     plot(Analog.Analog09.value)
                %                     plot(Analog.Analog10.value)
                %                     plot(Analog.Analog11.value)
                %                     plot(Analog.Analog12.value)
                %                     plot(Analog.Analog13.value)
                %                     plot(Analog.Analog14.value)
                %                     plot(Analog.Analog15.value)
                %                     plot(Analog.Analog16.value)
                %                     plot(Analog.Analog17.value)
                %                     plot(Analog.Analog18.value)
                %                     plot(Analog.Analog19.value)
                %                     plot(Analog.Analog20.value)
                %                     plot(Analog.Analog21.value)
                %                     plot(Analog.Analog22.value)
                %                     plot(Analog.Analog23.value)
                %                     plot(Analog.Analog24.value)
                %                     plot(Analog.Analog25.value)
                %                     plot(Analog.Analog26.value)
                %                     plot(Analog.Analog28.value)
                %                     plot(Analog.Analog29.value)
                %                     plot(Analog.Analog30.value)
                %
                %                 end

                % silence or auditory stimulation
                stF1= contains(filesList(fileInd(iFile)).name,'ilence');
                if stF1==1
                    silCondi=1;
                else
                    silCondi=0;
                end
                clear stF1


                %
                %                 if silCondi==0
                %                     %findpeaks(-tsout.Data,'MinPeakprominence',12,'MinPeakDistance',peakDAna,'Annotate','extents');
                [pks,locs,w,p] = findpeaks(Analog.Analog01.value,'MinPeakProminence',0.1,'MinPeakDistance',peakDAna);
                if mean(diff(locs))
                    fprintf('Warning!!! No analog signal! Tapping:  %s\n', subFolderNames{iPart})
                end
                %
                %                     diffLocs=diff(locs);
                %                     interRef=median(diffLocs);
                %                     iDeph=find(diffLocs>interRef+100 | diffLocs<interRef-100);
                %                     if iDeph(1)==1
                %                         iDeph=iDeph(2:end);
                %                     end
                %                     iAnaDeph=locs(iDeph+1);
                %
                %                     %Rank diff to classify trials
                %                     interMeas=diffLocs(iDeph);
                %                     DephMeas=interMeas-interRef;
                %                     interExp=[-interRef/4 interRef/4];
                %                     DephRank=nan(1,length(DephMeas));
                %                     for iii=1:length(DephMeas)
                %                         [~,DephRank(iii)] = min(reshape(abs(bsxfun(@minus,interExp,DephMeas(iii))),numel(interExp),[]));
                %                     end
                %
                %
                %                     % find reference deph =0
                %                     diffiDeph=diff(iDeph);
                %                     iiii=1;
                %                     if iDeph(1)>=10
                %                         iDephRef(iiii)=12;
                %                         iiii=iiii+1;
                %                     end
                %
                %                     for iii=1:length(diffiDeph)
                %                         if diffiDeph(iii)>=24
                %                             iDephRef(iiii)=iDeph(iii+1)-nbOfSteps;
                %                             iiii=iiii+1;
                %                         end
                %                     end
                %
                %                     if length(diffLocs)-iDeph(end)>=24
                %                         iDephRef(iiii)= iDeph(end)-nbOfSteps;
                %                     end
                %                     if exist('iDephRef') ==1
                %                         iAnaDephRef=locs(iDephRef);
                %                     end
                %
                %                     if exist('iAnaDephRef')==1
                %                         fprintf('No deph. (ref): %d, Deph. ref: %d\n', length(iAnaDephRef), length(DephRank));
                %                     else
                %                         fprintf('No deph. (ref): %d, Deph. ref: %d\n', 0, length(DephRank));
                %                     end
                %                     clear iDephRef

                %
                %
                %                     [iHSTOL,iHSTOR,iTOHSL,iTOHSR]=footKinCue(Data,1,11,2,12,Fs(1),5,0,0);
                %                     iHSTO=sortrows([iHSTOL;iHSTOR]);
                %                     iTOHS=sortrows([iTOHSL;iTOHSR]);
                %
                %
                %                     for iiDeph=1:2
                %                         iDephRank=find(DephRank==iiDeph);
                %                         if isempty(iDephRank)==0
                %                             for ii=1:length(iDephRank)
                %                                 % Find inex of steps close to deph detected
                %                                 % in anlog signal. The ration of sampling
                %                                 % frequency is considered in  bsxfun(@minus,iHSTO(:,1),iAnaDeph(iDephRank(ii))/2))
                %                                 % Sampling freq of analog =2222.2 Hz /
                %                                 % Sampling freq of insoles =1111.1 Hz
                %                                 [~,ind1] = min(reshape(abs(bsxfun(@minus,iHSTO(:,1),iAnaDeph(iDephRank(ii))/2)),numel(iHSTO(:,1)),[]));
                %                                 if abs(iHSTO(ind1)-iAnaDeph(iDephRank(ii))/2)<150
                %                                     l=length(iHSTO)-(ind1+5);
                %                                     if l>=0
                %                                         tmp1(ii,:)=iHSTO(ind1:ind1+5,3)';
                %                                     else
                %                                         tmp1(ii,:)=[iHSTO(ind1:ind1+5+l,3)',nan(1,abs(l))];
                %                                     end
                %                                     ind2=find(iTOHS(:,1)>iHSTO(ind1,1),1,'first');
                %                                     l=length(iTOHS)-(ind2+5);
                %                                     if l>=0
                %                                         tmp2(ii,:)=iTOHS(ind2:ind2+5,3)';
                %                                     else
                %                                         tmp2(ii,:)=[iTOHS(ind2:ind2+5+l,3)',nan(1,abs(l))];
                %                                     end
                %                                 end  %
                %                             end
                %                             if exist('tmp1')==1
                %                                 cuedGait(iPart,noCondi,iiDeph).CT=[cuedGait(iPart,noCondi,iiDeph).CT;tmp1];
                %                                 cuedGait(iPart,noCondi,iiDeph).CT(find(cuedGait(iPart,noCondi,iiDeph).CT(:,1)>0),:);
                %                                 clear tmp1
                %                             end
                %                             if exist('tmp2')==1
                %                                 cuedGait(iPart,noCondi,iiDeph).ST=[cuedGait(iPart,noCondi,iiDeph).CT;tmp2];
                %                                 cuedGait(iPart,noCondi,iiDeph).ST(find(cuedGait(iPart,noCondi,iiDeph).ST(:,1)>0),:);
                %                                 clear tmp2
                %                             end
                %                         end
                %
                %                     end
                %
                %                     % ref no deph
                %                     if exist('iAnaDephRef')==1
                %                         for ii=1:length(iAnaDephRef)
                %                             [~,ind1] = min(reshape(abs(bsxfun(@minus,iHSTO(:,1),iAnaDephRef(ii)/2)),numel(iHSTO(:,1)),[]));
                %                             if abs(iHSTO(ind1)-iAnaDephRef(ii)/2)<150
                %                                 l=length(iHSTO)-(ind1+5);
                %                                 if l>=0
                %                                     tmp1(ii,:)=iHSTO(ind1:ind1+5,3)';
                %                                 else
                %                                     tmp1(ii,:)=[iHSTO(ind1:ind1+5+l,3)',nan(1,abs(l))];
                %                                 end
                %                                 ind2=find(iTOHS(:,1)>iHSTO(ind1,1),1,'first');
                %                                 l=length(iTOHS)-(ind2+5);
                %                                 if l>=0
                %                                     tmp2(ii,:)=iTOHS(ind2:ind2+5,3)';
                %                                 else
                %                                     tmp2(ii,:)=[iTOHS(ind2:ind2+5+l,3)',nan(1,abs(l))];
                %                                 end
                %                             end
                %                         end
                %                         if exist('tmp1')==1
                %                             cuedGait(iPart,noCondi,3).CT=[cuedGait(iPart,noCondi,3).CT;tmp1];
                %                             cuedGait(iPart,noCondi,3).CT(find(cuedGait(iPart,noCondi,3).CT(:,1)>0),:);
                %                             clear tmp1
                %                         end
                %                         if exist('tmp2')==1
                %                             cuedGait(iPart,noCondi,3).ST=[cuedGait(iPart,noCondi,3).ST;tmp2];
                %                             cuedGait(iPart,noCondi,3).ST(find(cuedGait(iPart,noCondi,3).ST(:,1)>0),:);
                %                             clear tmp2 iAnaDephRef
                %                         end
                %                     end
                %
                %                 elseif silCondi==1
                %
                %                     % silCond if loop
                %
                %end
            end



        end



        %% Walking
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if noCondi==2

            peakDAna=90;%round(60/cad*250);


            for iFile=1:length(fileInd)
                %% Retrouver les dephasages du signal

                file2read=fullfile('../DATA/',subFolderNames{iPart},filesList(fileInd(iFile)).name);
                fprintf('Reading %s\n', file2read)


                if file2read(end)=='v'
                    % case CSV
                    Data = importCSV(file2read);
                else
                    % case .mat
                    load(file2read);
                end

                % silence or auditory stimulation
                stF1= contains(filesList(fileInd(iFile)).name,'ilence');
                if stF1==1
                    silCondi=1;
                else
                    silCondi=0;
                end
                clear stF1



                if silCondi==0
                    %                 if  isfield(Analog,'Channel01')==1
                    %                     tsin = timeseries(Analog.Channel01.value',1/Analog.Channel01.Rate:1/Analog.Channel01.Rate:length(Analog.Channel01.value)/Analog.Channel01.Rate);
                    %                 else
                    %                     tsin=timeseries(Analog.Analog01.value',1/Analog.Analog01.Rate:1/Analog.Analog01.Rate:length(Analog.Analog01.value)/Analog.Analog01.Rate);
                    %                 end
                    %                 tsout = detrend(tsin,'linear');
                    % tsout = detrend(Data(end,:),'linear');

                    %% Analog audio signal processing
                    ana=detrend(Data(end,:),'constant');
                    ana=downsample(ana,2);
                    %findpeaks(-tsout.Data,'MinPeakpromine nce',12,'MinPeakDistance',peakDAna,'Annotate','extents');
                    [pks,locs,w,p] = findpeaks(ana,'MinPeakProminence',0.03,'MinPeakDistance',peakDAna);


                    diffLocs=diff(locs);
                    interRef=median(diffLocs);
                    iDeph=find(diffLocs>interRef+100 | diffLocs<interRef-100);
                    if iDeph(1)==1
                        iDeph=iDeph(2:end);
                    end
                    iAnaDeph=locs(iDeph+1);

                    %Rank diff to classify trials
                    interMeas=diffLocs(iDeph);
                    DephMeas=interMeas-interRef;
                    interExp=[-interRef/4 interRef/4];
                    DephRank=nan(1,length(DephMeas));
                    for iii=1:length(DephMeas)
                        [~,DephRank(iii)] = min(reshape(abs(bsxfun(@minus,interExp,DephMeas(iii))),numel(interExp),[]));
                    end


                    % find reference deph =0
                    diffiDeph=diff(iDeph);
                    iiii=1;
                    if iDeph(1)>=20
                        iDephRef(iiii)=8;
                        iiii=iiii+1;
                    end

                    for iii=1:length(diffiDeph)
                        if diffiDeph(iii)>=12
                            iDephRef(iiii)=iDeph(iii+1)-6;
                            iiii=iiii+1;
                        end
                    end

                    if length(diffLocs)-iDeph(end)>=12
                        iDephRef(iiii)= iDeph(end)-6;
                    end
                    if exist('iDephRef','var') ==1
                        iAnaDephRef=locs(iDephRef);
                    end

                    if exist('iAnaDephRef','var')==1
                        fprintf('No deph. (ref): %d, Deph. ref: %d\n', length(iAnaDephRef), length(DephRank));
                    else
                        fprintf('No deph. (ref): %d, Deph. ref: %d\n', 0, length(DephRank));
                    end
                    clear iDephRef


                    %% Foot kinematics detection
                    [iHSTOL,iHSTOR,iTOHSL,iTOHSR]=footKinCue(Data,1,11,2,12,Fs(1),5,0,0);
                    iHSTO=sortrows([iHSTOL;iHSTOR]);
                    iTOHS=sortrows([iTOHSL;iTOHSR]);
                    % conversion of duration into ms
                    iHSTO(:,3)=iHSTO(:,3)*1000/Fs(1);
                    iTOHS(:,3)=iTOHS(:,3)*1000/Fs(1);

                    %% Pairing of foot kinematics with dephasing

                    for iiDeph=1:2
                        iDephRank=find(DephRank==iiDeph);
                        if isempty(iDephRank)==0
                            for ii=1:length(iDephRank)
                                % Find inex of steps close to deph detected
                                % in anlog signal. The ration of sampling
                                % frequency is considered in  bsxfun(@minus,iHSTO(:,1),iAnaDeph(iDephRank(ii))/2))
                                % Sampling freq of analog =2222.2 Hz /
                                % Sampling freq of insoles =1111.1 Hz

                                % find the closest step to the index of dephasing
                                [~,ind1] = min(reshape(abs(bsxfun(@minus,iHSTO(:,1),iAnaDeph(iDephRank(ii)))),numel(iHSTO(:,1)),[]));
                                % 150 threshold: number of samples between
                                % the dephased beat and its closest step
                                if abs(iHSTO(ind1)-iAnaDeph(iDephRank(ii)))<600
                                    l=length(iHSTO)-(ind1+nbOfSteps);
                                    if l>=0
                                        tmp1(ii,:)=iHSTO(ind1:ind1+nbOfSteps,3)';
                                        %relative phase
                                        iAnaDeph1stStep=iAnaDeph(iDephRank(ii));
                                        iAnaDephAllSteps=[iAnaDeph1stStep,locs(find(locs>iAnaDeph1stStep,nbOfSteps-1,'first'))];
                                        if length(iAnaDephAllSteps)==nbOfSteps
                                            for iii=1:nbOfSteps
                                                tmp3(ii,iii)=(iHSTO(ind1+iii-1,1)-iAnaDephAllSteps(iii))/interRef*2*pi;
                                            end
                                        end
                                        %                                     else
                                        %                                         tmp1(ii,:)=[iHSTO(ind1:ind1+nbOfSteps+l,3)',nan(1,abs(l))];
                                        %                                         %relative phase

                                        ind2=find(iTOHS(:,1)>iHSTO(ind1,1),1,'first');
                                        l=length(iTOHS)-(ind2+nbOfSteps);
                                        if l>=0
                                            tmp2(ii,:)=iTOHS(ind2:ind2+nbOfSteps,3)';
                                            %                                         else
                                            %                                             tmp2(ii,:)=[iTOHS(ind2:ind2+nbOfSteps+l,3)',nan(1,abs(l))];
                                        end
                                    end  %
                                end
                            end
                            if exist('tmp1','var')==1
                                cuedGait(iPart,noCondi,iiDeph).CT=[cuedGait(iPart,noCondi,iiDeph).CT;tmp1];
                                cuedGait(iPart,noCondi,iiDeph).CT(find(cuedGait(iPart,noCondi,iiDeph).CT(:,1)>0),:);
                                cuedGait(iPart,noCondi,iiDeph).relPh=[cuedGait(iPart,noCondi,iiDeph).relPh;tmp3];
                                clear tmp1
                            end
                            if exist('tmp3','var')==1
                                cuedGait(iPart,noCondi,iiDeph).relPh=[cuedGait(iPart,noCondi,iiDeph).relPh;tmp3];
                                clear tmp3
                            end
                            if exist('tmp2','var')==1
                                cuedGait(iPart,noCondi,iiDeph).ST=[cuedGait(iPart,noCondi,iiDeph).ST;tmp2];
                                cuedGait(iPart,noCondi,iiDeph).ST(find(cuedGait(iPart,noCondi,iiDeph).ST(:,1)>0),:);
                                clear tmp2
                            end
                        end

                    end

                    %% Ref no deph

                    if exist('iAnaDephRef','var')==1
                        for ii=1:length(iAnaDephRef)
                            % find the closest step to the index of dephasing
                            [~,ind1] = min(reshape(abs(bsxfun(@minus,iHSTO(:,1),iAnaDephRef(ii))),numel(iHSTO(:,1)),[]));
                            % 150 threshold: number of samples between
                            % the dephased beat and its closest step
                            if abs(iHSTO(ind1)-iAnaDephRef(ii))<600
                                l=length(iHSTO)-(ind1+nbOfSteps);
                                if l>=0
                                    tmp1(ii,:)=iHSTO(ind1:ind1+nbOfSteps,3)';
                                    %relative phase
                                    iAnaDeph1stStep=iAnaDephRef(ii);
                                    iAnaDephAllSteps=[iAnaDeph1stStep,locs(find(locs>iAnaDeph1stStep,nbOfSteps-1,'first'))];
                                    if length(iAnaDephAllSteps)==nbOfSteps
                                        for iii=1:nbOfSteps
                                            tmp3(ii,iii)=(iHSTO(ind1+iii-1,1)-iAnaDephAllSteps(iii))/interRef*2*pi;
                                        end
                                    end
                                    %                                 else
                                    %                                     tmp1(ii,:)=[iHSTO(ind1:ind1+nbOfSteps+l,3)',nan(1,abs(l))];
                                    %                                 end
                                    ind2=find(iTOHS(:,1)>iHSTO(ind1,1),1,'first');
                                    l=length(iTOHS)-(ind2+nbOfSteps);
                                    if l>=0
                                        tmp2(ii,:)=iTOHS(ind2:ind2+nbOfSteps,3)';
                                    else
                                        tmp2(ii,:)=[iTOHS(ind2:ind2+nbOfSteps+l,3)',nan(1,abs(l))];
                                    end
                                end
                            end
                        end
                        if exist('tmp1','var')==1
                            cuedGait(iPart,noCondi,3).CT=[cuedGait(iPart,noCondi,3).CT;tmp1];
                            cuedGait(iPart,noCondi,3).CT(find(cuedGait(iPart,noCondi,3).CT(:,1)>0),:);
                            clear tmp1
                        end
                        if exist('tmp3','var')==1
                            cuedGait(iPart,noCondi,3).relPh=[cuedGait(iPart,noCondi,3).relPh;tmp3];

                            clear tmp3
                        end

                        if exist('tmp2','var')==1
                            cuedGait(iPart,noCondi,3).ST=[cuedGait(iPart,noCondi,3).ST;tmp2];
                            cuedGait(iPart,noCondi,3).ST(find(cuedGait(iPart,noCondi,3).ST(:,1)>0),:);
                            clear tmp2 iAnaDephRef
                        end
                    end

                elseif silCondi==1

                    % silCond if loop

                end
            end



        end

        % end of condi loop
    end

    %% CLEANING and MEANS of OUTPUT VARIABLES

    for ii=1:3
        % Suppress NAN
        cuedGait(iPart,2,ii).ST=cuedGait(iPart,2,ii).ST(isnan(cuedGait(iPart,2,ii).ST(:,end))==0,:);
        cuedGait(iPart,2,ii).CT=cuedGait(iPart,2,ii).CT(isnan(cuedGait(iPart,2,ii).CT(:,end))==0,:);
        cuedGait(iPart,2,ii).relPh=cuedGait(iPart,2,ii).relPh(isnan(cuedGait(iPart,2,ii).relPh(:,end))==0,:);
        % Suppress Zeros
        cuedGait(iPart,2,ii).CT( all(~cuedGait(iPart,2,ii).CT,2), : ) = [];
        cuedGait(iPart,2,ii).ST( all(~cuedGait(iPart,2,ii).ST,2), : ) = [];
        cuedGait(iPart,2,ii).relPh( all(~cuedGait(iPart,2,ii).relPh,2), : ) = [];
        % Means
        cuedGait(iPart,2,ii).mCT=nanmean(cuedGait(iPart,2,ii).CT);
        cuedGait(iPart,2,ii).mST=nanmean(cuedGait(iPart,2,ii).ST);
        cuedGait(iPart,2,ii).mrelPh=circ_mean(cuedGait(iPart,2,ii).relPh);
    end

    % end of participant loop
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function donnesparticipants = read_part_Kq(workbookFile, sheetName, dataLines)

%% Input handling

% If no sheet is specified, read first sheet
if nargin == 1 || isempty(sheetName)
    sheetName = 1;
end

% If row start and end points are not specified, define defaults
if nargin <= 2
    dataLines = [4, 32];
end

%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 12);

% Specify sheet and range
opts.Sheet = sheetName;
opts.DataRange = "B" + dataLines(1, 1) + ":M" + dataLines(1, 2);

% Specify column names and types
opts.VariableNames = ["date", "nPart", "Code", "Nom", "Sexe", "Age", "BPMtapping", "FrquenceHz", "Priodem", "Cadence", "Poidskg", "Taillem"];
opts.VariableTypes = ["datetime", "double", "string", "string", "string", "double", "double", "double", "double", "double", "double", "double"];

% Specify variable properties
opts = setvaropts(opts, ["Code", "Nom", "Sexe"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Code", "Nom", "Sexe"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, "date", "InputFormat", "");

% Import the data
donnesparticipants = readtable(workbookFile, opts, "UseExcel", false);

for idx = 2:size(dataLines, 1)
    opts.DataRange = "B" + dataLines(idx, 1) + ":M" + dataLines(idx, 2);
    tb = readtable(workbookFile, opts, "UseExcel", false);
    donnesparticipants = [donnesparticipants; tb]; %#ok<AGROW>
end

end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function delsysDATA = importCSV(filename, dataLines)
%IMPORTFILE Import data from a text file
%  delsysDATA = IMPORTFILE(FILENAME) reads data from
%  text file FILENAME for the default selection.  Returns the data as a
%  table.
%
%  delsysDATA = IMPORTFILE(FILE, DATALINES) reads data
%  for the specified row interval(s) of text file FILENAME. Specify
%  DATALINES as a positive scalar integer or a N-by-2 array of positive
%  scalar integers for dis-contiguous row intervals.
%
%  Example:
%  delsysDATA = importfile("/Volumes/Data Hiva Oa/loicdamm/Library/CloudStorage/GoogleDrive-loic.damm@gmail.com/My Drive/Pro/cued gait 2.0 with code/DATA/003RIVA/delsys_003RIVA_shift2_20240117.csv", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 29-Jan-2024 13:19:11

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 42);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Xs", "FSRAdapter1FSRA1", "Xs1", "FSRAdapter1FSRB1", "Xs2", "FSRAdapter1FSRC1", "Xs3", "FSRAdapter1FSRD1", "Xs4", "FSRAdapter1ACCX1", "Xs5", "FSRAdapter1ACCY1", "Xs6", "FSRAdapter1ACCZ1", "Xs7", "FSRAdapter1GYROX1", "Xs8", "FSRAdapter1GYROY1", "Xs9", "FSRAdapter1GYROZ1", "Xs10", "FSRAdapter2FSRA2", "Xs11", "FSRAdapter2FSRB2", "Xs12", "FSRAdapter2FSRC2", "Xs13", "FSRAdapter2FSRD2", "Xs14", "FSRAdapter2ACCX2", "Xs15", "FSRAdapter2ACCY2", "Xs16", "FSRAdapter2ACCZ2", "Xs17", "FSRAdapter2GYROX2", "Xs18", "FSRAdapter2GYROY2", "Xs19", "FSRAdapter2GYROZ2", "Xs20", "AnalogInputAdapter3Analog3"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
delsysDATA = readtable(filename, opts);
delsysDATA=delsysDATA{:,2:2:42}';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [iHSTOL,iHSTOR,iTOHSL,iTOHSR]=footKinCue(FSR,channel_left_heel,channel_right_heel,channel_left_toe,channel_right_toe,fs,threshold,plot_yn,no_fig)
%[iHSTOL,iHSTOR,iTOHSL,iTOHSR]=footKinCue(iAnaDeph,calcL,toeL,calcR,toeR,6)
%[iHSL,iHSR,iTOL,iTOR,iHSTOL,iHSTOR,iTOHSL,iTOHSR]=footKinCue(FSR,channel_left_heel,channel_right_heel,channel_left_toe,channel_right_toe,fs,threshold,plot_yn,no_fig)

% FSR=Data;
% channel_left_heel=1;
% channel_right_heel=11;
% channel_left_toe=2;
% channel_right_toe=12;
% threshold=4;
% plot_yn=1;
% no_fig=1;
% no_subplot=1;
% fs=1111.11;

%%%%%%%%%%%%%%%%%%%%
%% Left
%%%%%%%%%%%%%%%%%%%%

P_heel=FSR(channel_left_heel,:);
P_toe=FSR(channel_left_toe,:);

f_cut=20; %  Cutoff frequency
step_length=200;
interval_ana=4000;

[b_butter, a_butter] = butter(5,2*f_cut/fs,'low');% f_uencies values are specified in normalized terms between 0.0 and 1.0,

ll=length(P_heel);
P_heel_filt=filtfilt(b_butter,a_butter,double(P_heel));
%P_toe_filt=filtfilt(b_butter,a_butter,double(P_toe));

%% detection of points near zero
inter_step=find(P_heel_filt<threshold);
transi=find(diff(inter_step)>step_length);


%% step detection

if isempty(transi)==0

    l=length(transi);
    iHSL=NaN(1,l-1);
    iTOL=NaN(1,l-1);

    for i=1:l
        if i==1
            tmp=round(median(inter_step(1):inter_step(transi(i))));
            tmp3=find(P_heel(tmp:tmp+interval_ana)>threshold,1,'first')+tmp;
            if isempty(tmp3)==1
                iHSL(i)=nan;
            else
                iHSL(i)=tmp3;
            end
        else
            tmp=round(median(inter_step(transi(i-1)+1):inter_step(transi(i))));
            if ll> tmp+interval_ana
                tmp3=find(P_heel(tmp:tmp+interval_ana)>threshold,1,'first')+tmp;
                if isempty(tmp3)==1
                    iHSL(i)=nan;
                else
                    iHSL(i)=tmp3;
                end
            end
        end

        if i==l
            tmp=round(median(inter_step(transi(i)+1):inter_step(end)));
            tmp3=-find(P_toe(tmp:-1:tmp-interval_ana)>threshold,1,'first')+tmp;
            if isempty(tmp3)==1
                iTOL(i)=nan;
            else
                iTOL(i)=tmp3;
            end
        else
            tmp=round(median(inter_step(transi(i)+1):inter_step(transi(i+1))));
            if tmp-interval_ana>=1
                tmp3=-find(P_toe(tmp:-1:tmp-interval_ana)>threshold,1,'first')+tmp;
                if isempty(tmp3)==1
                    iTOL(i)=nan;
                else
                    iTOL(i)=tmp3;
                end
            end
        end

    end

end

iHSL(isnan(iHSL))=[];
iTOL(isnan(iTOL))=[];
%
if plot_yn==1
    figure(no_fig)
    set(gcf,'color', 'white')
    subplot(2,1,1)
    plot(P_heel)
    hold on
    plot(P_toe)
    scatter(iHSL,P_heel(iHSL),'r')
    scatter(iTOL,P_toe(iTOL),'g')
end

%%%%%%%%%%%%%%%%%%%%
%% Right
%%%%%%%%%%%%%%%%%%%%
P_heel=FSR(channel_right_heel,:);
P_toe=FSR(channel_right_toe,:);

f_cut=20; %  Cutoff frequency
step_length=200;
interval_ana=4000;

[b_butter, a_butter] = butter(5,2*f_cut/fs,'low');% frequencies values are specified in normalized terms between 0.0 and 1.0,

ll=length(P_heel);
P_heel_filt=filtfilt(b_butter,a_butter,double(P_heel));

%% detection of points near zero
inter_step=find(P_heel_filt<threshold);
transi=find(diff(inter_step)>step_length);


%% step detection

if isempty(transi)==0

    l=length(transi);
    iHSR=NaN(1,l-1);
    iTOR=NaN(1,l-1);

    for i=1:l
        if i==1
            tmp=round(median(inter_step(1):inter_step(transi(i))));
            tmp3=find(P_heel(tmp:tmp+interval_ana)>threshold,1,'first')+tmp;
            if isempty(tmp3)==1
                iHSR(i)=nan;
            else
                iHSR(i)=tmp3;
            end
        else
            tmp=round(median(inter_step(transi(i-1)+1):inter_step(transi(i))));
            if ll> tmp+interval_ana
                tmp3=find(P_heel(tmp:tmp+interval_ana)>threshold,1,'first')+tmp;
                if isempty(tmp3)==1
                    iHSR(i)=nan;
                else
                    iHSR(i)=tmp3;
                end
            end
        end

        if i==l
            tmp=round(median(inter_step(transi(i)+1):inter_step(end)));
            tmp3=-find(P_toe(tmp:-1:tmp-interval_ana)>threshold,1,'first')+tmp;
            if isempty(tmp3)==1
                iTOR(i)=nan;
            else
                iTOR(i)=tmp3;
            end
        else
            tmp=round(median(inter_step(transi(i)+1):inter_step(transi(i+1))));
            if tmp-interval_ana>=1
                tmp3=-find(P_toe(tmp:-1:tmp-interval_ana)>threshold,1,'first')+tmp;
                if isempty(tmp3)==1
                    iTOR(i)=nan;
                else
                    iTOR(i)=tmp3;
                end
            end
        end

    end

end

iHSR(isnan(iHSR))=[];
iTOR(isnan(iTOR))=[];
%
if plot_yn==1
    figure(no_fig)
    set(gcf,'color', 'white')
    subplot(2,1,2)
    plot(P_heel)
    hold on
    plot(P_toe)
    scatter(iHSR,P_heel(iHSR),'r')
    scatter(iTOR,P_toe(iTOR),'g')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PAIRING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%[iHSTOL,iHSTOR,iTOHSL,iTOHSR]

%% PAIR TO & HS to identify flying phase
ii=1;
for i=1:length(iTOL)
    [~,iNearestiHS] = min(reshape(abs(bsxfun(@minus,iHSL,iTOL(i))),numel(iHSL),[]));
    %[~,iNearestiHS] = find(iHS>iTO(i),1,'first');

    if isempty(iNearestiHS)==0
        if iHSL(iNearestiHS)-iTOL(i)>10 && iHSL(iNearestiHS)-iTOL(i)<1000
            iTOHSL(ii,1:3)=[iTOL(i),iHSL(iNearestiHS),iHSL(iNearestiHS)-iTOL(i)];
            ii=ii+1;
        end
    end

end

ii=1;
for i=1:length(iTOR)
    [~,iNearestiHS] = min(reshape(abs(bsxfun(@minus,iHSR,iTOR(i))),numel(iHSR),[]));
    %[~,iNearestiHS] = find(iHS>iTO(i),1,'first');

    if isempty(iNearestiHS)==0
        if iHSR(iNearestiHS)-iTOR(i)>10 && iHSR(iNearestiHS)-iTOR(i)<1000
            iTOHSR(ii,1:3)=[iTOR(i),iHSR(iNearestiHS),iHSR(iNearestiHS)-iTOR(i)];
            ii=ii+1;
        end
    end

end

%% PAIR HS & TO to identify contact phase
ii=1;
for i=1:length(iHSL)
    %[~,iNearestiTO] = min(reshape(abs(bsxfun(@minus,iTO,iHS(i))),numel(iTO),[]));
    [~,iNearestiTO] = find(iTOL>iHSL(i),1,'first');

    if isempty(iNearestiTO)==0
        if iTOL(iNearestiTO)-iHSL(i)>10 && iTOL(iNearestiTO)-iHSL(i)<1000
            iHSTOL(ii,1:3)=[iHSL(i),iTOL(iNearestiTO),iTOL(iNearestiTO)-iHSL(i)];
            ii=ii+1;
        end
    end

end

ii=1;
for i=1:length(iHSR)
    %[~,iNearestiTO] = min(reshape(abs(bsxfun(@minus,iTO,iHS(i))),numel(iTO),[]));
    [~,iNearestiTO] = find(iTOR>iHSR(i),1,'first');

    if isempty(iNearestiTO)==0
        if iTOR(iNearestiTO)-iHSR(i)>10 && iTOR(iNearestiTO)-iHSR(i)<1000
            iHSTOR(ii,1:3)=[iHSR(i),iTOR(iNearestiTO),iTOR(iNearestiTO)-iHSR(i)];
            ii=ii+1;
        end
    end

end




end






