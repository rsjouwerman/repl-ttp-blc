%%% Baselinescript to extract baseline (2s prior to marker) and substract
%%% this from max after marker (CS-duration).
%%% For the ResDOPA (Mainz) data set
%%% 04/2019 Manuel Kuhn @ISN, Hamburg-Eppendorf
%%% adjust baseline corrections for multiverse baseline correction project
%%% 12/2020 Rachel Sjouwerman @ISN, Hamburg-Eppendorf

%%% 2021/03/16 I added postCSstart time 500 ms to correspond to the
%%% literature, as well as postCSend time -1500, 2000
%% Add path for pspm / spm
addpath('C:/Users/sjouwerman/Documents/Courses/PsPM_course/external/v4.3.0/')

%% Where to find data
presLogfiles = 'C:/Users/sjouwerman/Documents/repl-ttp-blc/data_raw/logfiles/';
path = 'C:/Users/sjouwerman/Documents/repl-ttp-blc/data_raw/matfiles/before_pspm_preprocessing/preprocessed/';
files = spm_select('FPList', path, 'tpspm_');

%% To be extracted data
% Allocate space
ROFLsubsData = nan(size(files,1),57);
allCSResponses = [];

%% Read in overview of specifics for different methods or loop over
for preCS = [1000, 2000, 3000, 4000, 9999] % 3000 ms is not actually in the literature
    for postCSstart = [0, 500, 1000, 2000, 3000] %
        for postCSend = [-1500, -1000, 0, 1000, 2000] % this actually stimulus duration plus this number
            for max1mean0 = [1, 0]
              
                
%                 sprintf('preCS %d ; postCS start %d ; post CS end 6000 - %d max (1) or mean (0) %d', ...
%                     preCS,postCSstart, postCSend, max1mean0);

                %% Study/Data specifics - please adjust
                samplingrate = 1000;
                % For rofl mininmal CS duration: 6000ms
                CSIntervall = 6000;
                
                postCSwindow_start = postCSstart;
                postCSwindow_end = CSIntervall + postCSend;
            
                    %postCSwindow = CSIntervall - postCSwindow_end;
                    % Define pre stimulus window for baseline
                    preCSwindow = preCS;
                    % Define postCS mean or max
                    postCSmax1mean0 = max1mean0; % for using max, else, use mean
                    
                    
                    
                    
                    
                    %% Study/Data specifics - please adjust
                    % samplingrate = 1000;
                    % % For rofl mininmal CS duration: 6000ms
                    % CSIntervall = 6000;
                    % %postCSwindow = CSIntervall;
                    % postCSwindow = 4000;
                    % % Define pre stimulus window for baseline
                    % preCSwindow = 999;
                    % % Define postCS mean or max
                    % postCSmax1mean0 = 1; % for using max, else, use mean
                    %% Get CS-duration infos here!
                    
                    
                    %% Extract data
                    % Extract Marker-wise for each subject and collect in Mat
                    for g = 1:size(files,1)
                        % Get subject info and load data
                        subFile = files(g,:);
                        subj = strsplit(subFile, '_');
                        subj = subj{end};
                        logFile = spm_select('FPList', presLogfiles, [subj(1:end-4) '.log']);
                        subj = subj(6:8);
                        subj = str2double(subj);
                        load(deblank(subFile));
                        
                        % Allocate to new variable and round marker time due to Samplingrate
                        scr = data{1}.data;
                        markerTimes = round(data{2}.data*samplingrate);
                        
                        %%% Run for each Marker
                        clear baseline TTPVal maxTTPVal
                        response = nan(1,57);
                        preCSval = nan(1,57);
                        postCSmean = nan(1,57);
                        postCSmax = nan(1,57);
                        postCSmax_onset = nan(1,57);
                        % Put subject number in first column
                        response(1) = subj;
                        preCSval(1) = subj;
                        postCSmean(1) = subj;
                        postCSmax(1) = subj;
                        postCSmax_onset(1) = subj;
                        % Fill the following columns with the respective response
                        for m = 1:length(markerTimes)
                            % Calculate Baseline
                            baseline = mean(scr((markerTimes(m)-preCSwindow):markerTimes(m)));
                            % Get SCR during CS-duration intervall after marker
                            TTPVal = scr((markerTimes(m)+postCSwindow_start):(markerTimes(m)+postCSwindow_end));
                            % Find max for this intervall and corresponding timepoint
                            [maxTTPVal, onset_maxTTPval] = max(TTPVal);
                            % Find mean for this intervall
                            meanTTPVal = mean(TTPVal);
                            % Define whether max or mean should be substracted
                            if true(postCSmax1mean0)
                                postCSval = maxTTPVal;
                            else
                                postCSval = meanTTPVal;
                            end
                            
                            % are there multiple peaks within the postCS time window
                            [pks, locs] = findpeaks(TTPVal);
                            
                            % subtract baseline from max or mean of intervall
                            response(m+1) = postCSval - baseline;
                            % extract other variables of interest: 1. preCSval, i.e., baseline
                            preCSval(m+1) = baseline;
                            % extract other variables of interest: 1. postCSval, i.e., mean/max
                            postCSmean(m+1) = meanTTPVal;
                            % extract other variables of interest: 1. postCSval, i.e., mean/max
                            postCSmax(m+1) =  maxTTPVal;
                            % extract other variables of interest: 1. onset postCSval
                            postCSmax_onset(m+1) = onset_maxTTPval;
                            
                            
                            %figure; plot(scr((markerTimes(m)-1999):(markerTimes(m)+CSIntervall)))
                        end
                        
                        % Collect subject's responses in its respective row
                        ROFLsubsData(g,:) = response;
                        
                        
                        %%% Add here CSP/CSM/US Conditions
                        
                        filename = logFile;
                        delimiter = '\t';
                        startRow = 4;
                        formatSpec = '%s%s%s%s%s%s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
                        fileID = fopen(filename,'r');
                        textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false, 'EndOfLine', '\r\n');
                        dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'ReturnOnError', false);
                        fclose(fileID);
                        eventData = [dataArray{1:end-1}];
                        clearvars filename delimiter startRow formatSpec fileID dataArray ans;
                        
                        
                        
                        %%%  Remove VAS Face Labels to avoid confusion in CS Onsets
                        vas1Ind = strfind(eventData(:,4), 'vas');
                        vas1Ind = find(not(cellfun('isempty',vas1Ind)));
                        eventData(vas1Ind,4) = 'VAS_noLabel';
                        
                        %%% Get Indices and Timings for Shocks and NoShocks
                        shockInd = strfind(eventData(:,4), 'Shock');
                        shockInd = find(not(cellfun('isempty',shockInd)));
                        shockInd = shockInd(1:3:end);
                        shockInd(:,2) = 999;
                        shockInd(:,3) = 1:length(shockInd);
                        
                        %%% Get Indices and Timing of CSs (Faces)
                        face1Ind = strfind(eventData(:,4), 'Face1');
                        face1Ind = find(not(cellfun('isempty',face1Ind)));
                        face1Ind(:,2) = 1;
                        face1Ind(:,3) = 1:length(face1Ind);
                        
                        face2Ind = strfind(eventData(:,4), 'Face2');
                        face2Ind = find(not(cellfun('isempty',face2Ind)));
                        face2Ind(:,2) = 2;
                        face2Ind(:,3) = 1:length(face2Ind);
                        
                        mFSEventLine = sortrows([face1Ind; face2Ind; shockInd ],1);
                        mFSFace1Ind = find(mFSEventLine(:,2) == 1);
                        mFSFace2Ind = find(mFSEventLine(:,2) == 2);
                        mFSShockInd = find(mFSEventLine(:,2) == 999);
                        
                        %%%% Which one is CSP
                        if sum(ismember(mFSShockInd-1,mFSFace1Ind(8:21))) == 14
                            cspInd = find((mFSEventLine(:,2) == 1));
                            csmInd = find((mFSEventLine(:,2) == 2));
                        elseif sum(ismember(mFSShockInd-1,mFSFace2Ind(8:21))) == 14
                            cspInd = find((mFSEventLine(:,2) == 2));
                            csmInd = find((mFSEventLine(:,2) == 1));
                        else
                            error('no correct CSP/CSM assignment');
                        end
                        
                        
                        mFSEventLine(:,4) = NaN;
                        mFSEventLine(find(mFSEventLine(:,2) == 999),4) = 999;
                        mFSEventLine(cspInd,4) = 1;
                        mFSEventLine(csmInd,4) = -1;
                        mFSEventLine(:,5) = response(2:end);
                        %add other values of interest
                        mFSEventLine(:,6) = preCSval(2:end);
                        mFSEventLine(:,7) = postCSmean(2:end);
                        mFSEventLine(:,8) = postCSmax(2:end);
                        mFSEventLine(:,9) = postCSmax_onset(2:end);
                        
                        
                        %subCSResponses = [repmat(subj,size(mFSEventLine,1),1) mFSEventLine(:,3:end)];
                        
                        subCSResponses = [repmat(subj,size(mFSEventLine,1),1) ...
                            repmat(postCSwindow_start,size(mFSEventLine,1),1) ...
                            repmat(postCSwindow_end,size(mFSEventLine,1),1) ...
                            repmat(preCSwindow,size(mFSEventLine,1),1) ...
                            repmat(postCSmax1mean0,size(mFSEventLine,1),1) ...
                            mFSEventLine(:,3:end)];
                        
                        allCSResponses = [allCSResponses; subCSResponses];
                        
                        
                    end
           
                
                
            end % these ar the ends of the pre post (start and end) time interval windows
        end % and the mean or max
    end
end
figure;
plot(mean(ROFLsubsData(:,2:end),1))



%% Save data
save('C:/Users/sjouwerman/Documents/repl-ttp-blc/data/baseline_corrected_multiverse_all_subjects.mat', 'allCSResponses');
csvwrite('C:/Users/sjouwerman/Documents/repl-ttp-blc/data/baseline_corrected_multiverse_all_subjects.csv', allCSResponses);

%%save('/home/kuhn/projects/pspm_Vergleich/Day1_ROFL_baselineTTP.mat', 'allCSResponses');
%%csvwrite('/home/kuhn/projects/pspm_Vergleich/Day1_ROFL_baselineTTP.csv', allCSResponses);
