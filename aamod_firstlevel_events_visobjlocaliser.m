% Generate names/onsets/durations struct array and save (where?). You must
% provide the correct sessions, e.g. by branching. And all sessions must
% include 'visobjlocaliser' to help us find the correct session number
function [aap,resp]=aamod_firstlevel_events_visobjlocaliser(aap,task,subj)

resp = '';

switch task
    case 'report'

    case 'doit'
        % find subject name
        subname = aap.acq_details.subjects(subj).mriname;
        % find exp dir
        expdir = fileparts(which('exp_visobjlocaliser'));
        infile = fullfile(expdir,'subjects',subname,...
            'data_exp_visobjlocaliser','subdata.mat');
        % load events
        subdata = loadbetter(infile);
        data = subdata2aa(subdata,aap,...
            aap.tasklist.currenttask.settings.sessiontarget,1);
        nsess = length(aap.acq_details.selected_sessions);
        assert(nsess == length(data),'subdata / session mismatch');
        spmpath = aas_getfiles_bystream(aap,subj,'firstlevel_spm');
        load(spmpath);
        for sess = 1:nsess
            % parse into nod
            NOD = struct('names',...
                {{'faces','places','objects','scrambled','responses'}},...
                'onsets',{repmat({[]},[1,5])},'durations',...
                {repmat({[]},[1,5])});
            sessres = data(sess).res;
            % we had 36 presentations per block so these are the block
            % onsets
            blocks = sessres(1:36:end,:);
            for b = 1:length(blocks)
                bname = blocks{b,1};
                if strfind(bname,'scrambled');
                    category = 'scrambled';
                % this is how bad Matlab is at handling strings
                elseif any(cell2mat(cellfun(@strfind,...
                        repmat({bname},[1 4]),...
                        {'Indoor','Place','House','Scene'},...
                        'uniformoutput',false)))
                    category = 'places';
                elseif ~isempty(str2num(bname(1)))
                    category = 'faces';
                else
                    % this is a bit unsafe but there's just too many
                    % objects and I cannae be bothered
                    category = 'objects';
                end
                % find the category
                catind = strfindcell(NOD.names,category);
                NOD.onsets{catind} = [NOD.onsets{catind} ...
                    blocks{b,5}];
                NOD.durations{catind} = [NOD.durations{catind} 16];
            end
            % add responses
            respvec = logical(cell2mat(sessres(:,3)));
            respons = cell2mat(sessres(respvec,5))';
            NOD.onsets{5} = respons;
            % impulses
            NOD.durations{5} = 0;
            % add to spm
            for c = 1:length(NOD.names)
                SPM.Sess(sess).U(c) = struct(...
                    'ons',NOD.onsets{c},... 
                    'dur',NOD.durations{c},...
                    'name',{NOD.names(c)},...
                    'P',struct('name','none'));
            end
            % Assume no nuisance regressors (but if you put some in
            % already that should work too)
            if ~isfield(SPM.Sess(sess),'C') || ~isstruct(SPM.Sess(sess).C)
                SPM.Sess(sess).C.C = [];
                SPM.Sess(sess).C.name = {};
            end
        end
        % save and describe
        save(spmpath,'SPM');
        aap = aas_desc_outputs(aap,subj,'firstlevel_spm',spmpath);
    case 'checkrequirements'
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
