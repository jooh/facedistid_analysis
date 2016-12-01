% Generate names/onsets/durations struct array and save. You must
% provide the correct sessions, e.g. by branching. 
% Now also upcasts the SPM to look like 16 sessions (so each sub-run is
% treated as a different run)
function [aap,resp]=aamod_firstlevel_events_oneback(aap,task,subj)

resp = '';

switch task
    case 'report'

    case 'doit'
        pptpath = false;
        if isempty(which('KbCheck'))
            oldpath = path;
            start_psychtoolbox;
            pptpath = true;
        end

        % find subject name
        subname = aap.acq_details.subjects(subj).mriname;
        % find exp dir
        expdir = fileparts(which('exp_oneback'));
        infile = fullfile(expdir,'subjects',subname,...
            'data_exp_oneback','subdata.mat');
        % load events
        subdata = loadbetter(infile);
        ts = aap.tasklist.currenttask.settings;
        % 2 data entries per session (16)
        data = subdata2aa(subdata,aap,ts.sessiontarget,2);
        % should be 8
        nsess = length(aap.acq_details.selected_sessions);
        % should be 16
        ndata = length(data);
        assert(nsess == (ndata/2),'subdata / session mismatch');
        spmpath = aas_getfiles_bystream(aap,subj,'firstlevel_spm');
        load(spmpath);
        % get time0 (same for 2 subruns)
        starttime = NaN([1 ndata]);
        for odds = 1:2:length(data)
            starttime(odds:odds+1) = data(odds).res.timestart;
        end
        % save the old SPM for indexing
        SPMorg = SPM;
        tr = SPM.xY.RT;
        % define converters for going become time in s and volumes.
        % (nb, time0 = volume 1. Since time is 0-indexed by
        % definition we must add 1 here to get to 1-indexing.)
        time2tr = @(time) 1 + (time/tr);
        tr2time = @(trtime) (trtime-1) * tr;

        % store time offset for events (0 if trim first, otherwise not)
        % (used to compensate timings if we have taken away the scans from
        % first half of run)
        offtime = zeros([1 ndata]);

        % if we aren't splitting, combine the subruns to a single data per
        % run (timings should all be relative to odds anyway
        if ~ts.splitsubruns
            newdata = data;
            for d=2:2:ndata
                % concatenate all the relevant fields
                newdata(d/2).res.trials = ...
                    [data(d-1).res.trials data(d).res.trials];
                for c = 1:length(data(1).res.conditions)
                    newdata(d/2).res.conditions(c).result = ...
                        [data(d-1).res.conditions(c).result;
                        data(d).res.conditions(c).result];
                end
            end
            % strip off extra datas
            starttime = starttime(1:2:ndata);
            ndata = ndata/2;
            data = newdata(1:ndata);
        else
            % cell array of path names
            volnames = cellstr(SPM.xY.P);
            % clear it out, to be reconstructed...
            SPM.xY.P = [];
            % amount of padding to apply 
            % (now less to get a better baseline estimate)
            padvol = 0; %2/tr;
        end
        % indices into the original SPM
        sessinds = repinds(1,nsess,1+ts.splitsubruns);

        for d = 1:ndata
            % find subdata 
            datres = data(d).res;
            % lots of acrobatics to split each SPM.mat run in 2
            if ts.splitsubruns
                % indices to volumes belonging to this run
                volinds = strfindcell(volnames, ...
                    aap.acq_details.sessions(...
                    aap.acq_details.selected_sessions(sessinds(d))).name);
                % keep track of where this adventure started
                vol1 = volinds(1);
                % need to offset time 0, only consider last n scans
                if iseven(d)
                    % reduce padding by half to get the same number of volumes
                    % for first and second splits
                    padding = padvol/2;
                    % need to make up a new session
                    % find the time of the first null trial
                    firsttime = datres.trials(1).time(1)-starttime(d);
                    % convert to tr
                    firstvol = time2tr(firsttime);
                    % add padding (null trial is 5 s so this should put first
                    % event at around 5 s)
                    firstvolind = floor(firstvol+padding);
                    % now positively identify wanted inds instead of removing
                    volinds = volinds(firstvolind:end);
                    % this is now time 0
                    offtime(d) = tr2time(firstvolind);
                else
                    padding = padvol;
                end
                
                % regardless of whether we have trimmed the beginning or not,
                % trim the end
                % AH HA! this is bad because endtime is inaccurate.
                % instead, use soa field to manually calculate somewhat
                % less inaccurate actual endtime
                lasttime = (datres.trials(end).time(1) + ...
                    datres.trials(end).condition.soa) - ...
                    (starttime(d) + offtime(d));
                % convert to volumes (nb floating point)
                lastvol = time2tr(lasttime);
                % find the last volume
                lastvolind = floor(lastvol-padding);
                volinds = volinds(1:lastvolind);
                % set nscan and add the volumes to the reconstructed xY.P
                SPM.nscan(d) = length(volinds);
                SPM.xY.P = [SPM.xY.P; SPMorg.xY.P(volinds,:)];
                % now we can parse the events for this run
                NOD = parseconditions(datres.conditions(1:end-1),...
                    datres.timeind,(starttime(d)+offtime(d)),...
                    ts.collapseview);
            else
                % nice and easy
                NOD = parseconditions(datres.conditions(1:end-1),...
                    datres.timeind,starttime(d),ts.collapseview)
                % add feedback event (from end of second null trial to
                % start of third)
                feedstart = datres.conditions(end).result(2).time(1)+ ...
                    datres.conditions(end).soa;
                % should be just over 5 s
                feeddur = datres.conditions(end).result(3).time(1) - ...
                    feedstart;
                feedons = feedstart - starttime(d);
                NOD.names{end+1} = 'feedback';
                NOD.onsets{end+1} = feedons;
                NOD.durations{end+1} = feeddur;
            end
            % add responses
            NOD.names{end+1} = 'responses';
            NOD.onsets{end+1} = cell2mat([datres.trials.responsetime]) ...
                - (starttime(d) + offtime(d));
            NOD.durations{end+1} = 0;
            % add to SPM
            for c = 1:length(NOD.names)
                % check that nothing is outside the available scans for the
                % sub-sampled runs
                assert(all(NOD.onsets{c}>0 & ...
                    NOD.onsets{c}<tr2time(SPM.nscan(d))),...
                    'onsets outside of run duration');
                SPM.Sess(d).U(c) = struct(...
                    'ons',NOD.onsets{c},... 
                    'dur',NOD.durations{c},...
                    'name',{NOD.names(c)},...
                    'P',struct('name','none'));
            end
            % Assume no nuisance regressors (but if you put some in
            % already that should work too)
            if isfield(SPMorg,'Sess') && isfield(...
                    SPMorg.Sess(sessinds(d)),'C')
                % find the run indices by subtracting the first vol in the
                % run (this may fail if you use dummy scans)
                runinds = volinds - (vol1-1);
                % also parse the covariates
                C = SPMorg.Sess(sessinds(d)).C.C(runinds,:);
                SPM.Sess(d).C = struct('C',C,'name',...
                    {SPMorg.Sess(sessinds(d)).C.name});
                % strip regressors that ended up empty after session split
                badinds = all(SPM.Sess(d).C.C==0,1);
                SPM.Sess(d).C.C(:,badinds) = [];
                SPM.Sess(d).C.name(badinds) = [];
            else
                SPM.Sess(d).C.C = [];
                SPM.Sess(d).C.name = {};
            end
        end % / d in ndata
        % need to make sure the high pass filter has been set correctly
        SPM.xX.K = repmat(SPM.xX.K(1),[1 ndata]);

        assert(length(unique(cellstr(SPM.xY.P)))==size(SPM.xY.P,1),...
            'duplicate volumes detected!');

        save(spmpath,'SPM');
        aap = aas_desc_outputs(aap,subj,'firstlevel_spm',spmpath);

        if pptpath
            % remove psychtoolbox again since it messes up printing
            path(oldpath);
        end
    case 'checkrequirements'
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end

function NOD = parseconditions(conditions,timeind,normtime,collapseview)

ncon = length(conditions);
NOD = struct('names',{{conditions.name}},...
    'onsets',{repmat({[]},[1,ncon])},...
    'durations',{repmat({[]},[1,ncon])});
getonsets = @(x) x(timeind) - normtime;
for c = 1:ncon
    NOD.onsets{c} = cellfun(getonsets,{conditions(c).result.time});
    % now model duration correctly
    NOD.durations{c} = repmat(2,size(NOD.onsets{c}));
end
if collapseview
    for c = 1:12
        [NOD.onsets{c},sortind] = sort([NOD.onsets{c} NOD.onsets{c+12}]);
        tdur = [NOD.durations{c} NOD.durations{c+12}];
        NOD.durations{c} = tdur(sortind);
    end
    NOD.names(13:end) = [];
    NOD.durations(13:end) = [];
    NOD.onsets(13:end) = [];
end
