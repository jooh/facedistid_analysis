% Create predictor RDMs for idanddistinct study.
% [aap,resp]=aamod_pilab_rsapredictors_oneback(aap,task,subj)
function [aap,resp]=aamod_pilab_rsapredictors_oneback(aap,task,subj)

resp = '';

switch task
    case 'report'

    case 'doit'
        % find subject name
        subname = aap.acq_details.subjects(subj).mriname;
        % find exp dir
        expdir = fileparts(which('exp_oneback'));
        subexp = fullfile(expdir,'subjects',subname);

        % process stimulusspace - basic predictions across full matrix
        stimpath = fullfile(subexp,'nose_stim_spherevectors.mat');
        ss = loadbetter(stimpath);
        rdms = struct('name','coef_eucdist','RDM',...
            ss.rdmbyattribute('shapecoef','euclidean'));
        rdms(2) = struct('name','coef_angdist','RDM',...
            ss.rdmbyattribute('shapecoef','angdist'));
        rdms(3) = struct('name','coef_eccdist','RDM',...
            ss.rdmbyattribute('shapecoef','eccdist'));
        rdms(4) = struct('name','vid_corrdist','RDM',...
            ss.rdmbyattribute('videoview','corr'));
        % indices corresponding to basic predictors
        baseinds = 1:4;
        clear ss

        % then consider different sub-spaces:

        % view-specific variants
        diagonal = logical(eye(24));
        blockf = false(12,12);
        blockt = true(12,12);
        % parts to mask for right-specific analysis
        notright = [blockf blockt; blockt blockt];
        notright(diagonal) = false;
        viewvariant.name = 'right';
        viewvariant.nanmask = notright;
        % parts to mask for left-specific
        notleft = [blockt blockt; blockt blockf];
        notleft(diagonal) = false;
        viewvariant(2).name = 'left';
        viewvariant(2).nanmask = notleft;
        % across view
        notacross = [blockt, blockf; blockf blockt];
        notacross(diagonal) = false;
        viewvariant(3).name = 'across';
        viewvariant(3).nanmask = notacross;

        for r = baseinds
            for v = 1:length(viewvariant)
                rd = rdms(r).RDM;
                rd(viewvariant(v).nanmask) = NaN;
                rdms(end+1) = struct('name',sprintf('view_%s_%s',...
                    viewvariant(v).name,rdms(r).name),'RDM',rd);
            end
        end

        % distinctiveness-specific models
        subs = repmat([true(1,3); true false(1,2); true false(1,2)],[8 8]);
        typicals = circshift(subs,[1 1]);
        caricatures = circshift(typicals,[1 1]);
        dist.subs = ~typicals & ~caricatures;
        dist.typicals = ~subs & ~caricatures;
        dist.caricatures = ~subs & ~typicals;
        dist.notsubs = ~subs;
        for rdind = 1:length(rdms)
            % distspecific analysis only really makes sense for angdist
            if isempty(strfind(rdms(rdind).name,'angdist'))
                continue
            end
            for fn = fieldnames(dist)'
                rdang = rdms(rdind).RDM;
                rdang(~dist.(fn{1})) = NaN;
                % make sure we don't nan out the diagonal since this can
                % lead to unpredictable pilab behaviour (mats being
                % interpreted as stacked vecs)
                rdang(diagonal) = 0;
                rdms(end+1) = struct('name',sprintf(...
                    'dspec_%s_%s',fn{1},rdms(rdind).name),...
                    'RDM',rdang);
            end
        end
        nrelevant = length(rdms);

        pairpath = fullfile(subexp,'data_exp_facepairs','subdata.mat');
        subdata = loadbetter(pairpath);
        assert(length(subdata)==16 || strcmp(lower(subname),'be'),...
            'incorrect number of facepairs sessions')
        rdsum = zeros(12,12);
        for sess = 1:length(subdata)
            rdsum = rdsum + subdata(sess).res.rdm;
        end
        rdms(end+1) = struct('name','facepairs','RDM',rdsum);
        clear subdata
        
        % score questionnaire data
        ratepath = fullfile(subexp,'data_exp_vidrate','subdata.mat');
        subdata = loadbetter(ratepath);
        % get items
        items = readquestcsv(fullfile(expdir,'exp_imagerate_qitems.csv'));
        % one catch - apparently LT only had a single session of the task
        % (everyone else had 2 I think)
        for it = 1:length(items)
            rating = zeros(12,1);
            for sess = 1:length(subdata)
                keys = subdata(sess).res.validkeys;
                score = subdata(sess).res.score;
                [junk,sessrating] = ismember(score(it+1).respk,keys);
                rating = rating + sessrating';
            end
            meanrating = rating / length(subdata);
            % extract construct name
            itname = textscan(items(it).label_high,'%s');
            itname = sprintf('vidrate-%s',itname{1}{end});
            rdms(end+1) = struct('name',itname,'RDM',squareform(pdist(...
                meanrating,'meandist')));
        end
        % Also do patterns across constructs - basic emotions and
        % dominance/trustworthiness
        compounds = struct('name',{'basicemotion','todorov'},...
            'inds',{5:10,11:12});
        for c = 1:length(compounds)
            ninds = length(compounds(c).inds);
            rating = zeros(12,ninds);
            for it = 1:ninds
                compind = compounds(c).inds(it);
                for sess = 1:length(subdata)
                    keys = subdata(sess).res.validkeys;
                    score = subdata(sess).res.score;
                    [junk,sessrating] = ismember(score(compind+1).respk,...
                        keys);
                    rating(:,it) = rating(:,it) + sessrating';
                end
            end
            meanrating = rating / length(subdata);
            rdms(end+1) = struct('name',sprintf('vidrate-%s',...
                compounds(c).name),'RDM',squareform(pdist(...
                meanrating,'euclidean')));
        end
        % upcast if necessary
        for n = 1:length(rdms)
            if size(rdms(n).RDM,1)==12
                rdms(n).RDM = repmat(rdms(n).RDM,[2 2]);
            end
        end
        % binary viewpoint RDM
        rdms(end+1) = struct('name','view','RDM',squareform(pdist(...
            [ones(12,1); zeros(12,1)],'meandist')));
        % Later: make behavioural / ET RDMs too
        % DEBUG: only key RDMs for now
        rdms = rdms([length(rdms) 1:nrelevant]);
        % save and describe
        pidir = fullfile(aas_getsubjpath(aap,subj),'pilab');
        mkdirifneeded(pidir);
        outpath = fullfile(pidir,'idanddistinct_rsapredictors.mat');
        save(outpath,'rdms');
        figdir = fullfile(pidir,'figures');
        mkdirifneeded(figdir);
        aap = aas_desc_outputs(aap,subj,'pilab_rsapredictors',outpath);
        % Also make quick diagnostic figures
        % get stimuli
        spath = aas_getfiles_bystream(aap,subj,'pilab_stimuli');
        stimuli = loadbetter(spath);
        for r = 1:length(rdms)
            % base rdm
            F = plotrdms(rdms(r),'labels',{stimuli.image},'nrows',3,...
                'titles',rdms(r).name);
            printstandard(fullfile(figdir,['predictor_rdm_' rdms(r).name]));
            close(F);
            % mds - wait for non-broken imaging system
            %%if ~any(isnan(rdms(r).RDM(:)))
                %F = plotrdms(rdms(r),'labels',{stimuli.image},'domds',true,...
                    %'dordm',false,'titles',rdms(r).name,'imagealpha',...
                    %{stimuli.alpha},'alphacolor',255);
                %printstandard(fullfile(figdir,['predictor_mds_' rdms(r).name]));
                %close(F);
            %else
                %fprintf('skipping mds for %s due to nans in matrix\n',...
                    %rdms(r).name);
            %end
            % image pixels
            im = intensity2rgb(rdms(r).RDM,jet(1e3));
            imwrite(im,fullfile(figdir,sprintf('predictor_rawrdm_%s.png',...
                rdms(r).name)),'PNG');
        end
        % show second-order RSM
        F = figure;
        % as big as it will go...
        set(F,'units','pixels','position',get(0,'screensize'));
        rsm = corr(asrdmvec(rdms),'rows','pairwise','type','spearman');
        F = plotrdms(rsm,'labels',stripbadcharacters({rdms.name},' '),...
            'docb',true,'colorbarargs',{'label','spearman rho'},...
            'titles','second-order similarity','fighand',F,...
            'rotatelabels',90,'cmap',jet(1e3));
        printstandard(fullfile(figdir,'predictor_so_rsm'));
        close(F);
        % mds - nah, too many nans
        %F = plotrdms(1-rsm,'labels',{rdms.name},'domds',true,'titles',...
            %'second-order similarity','dordm',false);
        %printstandard(fullfile(figdir,'predictor_so_mds'));
        %close(F);
    case 'checkrequirements'
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
