% evaluate effects of running inference on foldwise performance, even though the
% NFold CV scheme introduces a dependence (the test on one fold is part of the
% model in all the other folds).
%
% We test by generating an additional validation sample and testing each
% foldwise model twice: once on the (dependent) test set within the main sample,
% and again on an (independent) test set of equal size from the validation
% sample.
% 
% The simulation has a number of fixed parameters at the top of the function,
% which have been set to match the real data as close as possible.
%
% This code was developed and tested on Matlab R2013a.
%
% % INPUTS:
% recompute (false): rerun simulation even if existing result exists
% nsim (1000): number of iterations
% debug (false): end function in debugger
%
% 2016-12-08 J Carlin
%
% test_cvdependence(recompute,nsim,debug)
function test_cvdependence(recompute,nsim,debug)

% input check
if ~exist('recompute','var') || isempty(recompute)
    recompute = false;
end
% check if we can re-generate plots from an existing result
fundir = fileparts(mfilename('fullpath'));
cache = fullfile(fundir,'test_cvdependence_cache.mat');
hascache = exist(cache,'file');
if ~recompute && ~hascache
    fprintf('no previous result exists, recomputing.\n');
    recompute = true;
end
if ~exist('nsim','var') || isempty(nsim)
    nsim = 1000;
end
if ~exist('debug','var') || isempty(debug)
    debug = false;
end

% fixed parameters
nstim = 24;
nsub = 10;
nfeat = 200;
nreg = 4;
siglevel = linspace(0,.3,6);
nsiglevel = numel(siglevel);

if recompute
    % signal level for model 1, signal level for model 2, 2 models + their
    % difference
    outdim = [nsiglevel,nsiglevel,3];
    means_rcv = cell(nsim,1);
    means_rvalid = cell(nsim,1);
    stdev_rcv = cell(nsim,1);
    stdev_rvalid = cell(nsim,1);
    fprintf('starting simulation (%d iterations)...\n',nsim);
    tstart = tic;
    % awkwardness courtesy of matlab parfor limitations - we assign to a cell
    % array inside the parfor and then reconstruct the full 4D result after it's
    % done.
    parfor simind = 1:nsim
        % two models - now we recook on each parfor iteration
        model = struct;
        for m = 1:2
            % arbitrary model - each model feature has a positive loading on each stim
            model(m).regs = rand(nstim,nreg,'single');
            % model RDMs - squared euclidean (abs) distances in each
            model(m).modelrdv = arrayfun(@(x)pdist(model(m).regs(:,x),...
                @(a,b)abs(bsxfun(@minus,a,b)).^2)',1:nreg,'uniformoutput',0);
            model(m).modelrdv = horzcat(model(m).modelrdv{:});
            % arbitrary mapping to features - each feature is some linear combination of
            % regs
            model(m).weights = rand(nreg,nfeat,'single');
            % scale signal to mean=0, st dev 1 over features to match noise
            model(m).signal = zscore(model(m).regs * model(m).weights,[],2);
        end
        % cook noise - separate for each subject (and for the
        % independent validation sample)
        noise = randn(nstim,nfeat,nsub*2,'single');
        % awkwardness courtesy of matlab parfor limitations
        sim_means_rcv = NaN(outdim);
        sim_means_rvalid = NaN(outdim);
        sim_stdev_rcv = NaN(outdim);
        sim_stdev_rvalid = NaN(outdim);
        % run on various SNR combinations for the two models
        for sigind1 = 1:nsiglevel
            signal1 = model(1).signal * siglevel(sigind1);
            for sigind2 = 1:nsiglevel
                signal2 = model(2).signal * siglevel(sigind2);
                % bsxfun to expand third dim (independent noise for each subject)
                data = bsxfun(@plus,signal1+signal2,noise);
                % build euclidean.^2 RDV for each subject's data (stacked in second
                % dim)
                datardv = arrayfun(@(x)pdist(data(:,:,x))'.^2,1:nsub*2,...
                    'uniformoutput',0);
                datardv = horzcat(datardv{:});
                % fit multi reg RSA model to each subject and model
                rcv = NaN([nsub,2]);
                rvalid = rcv;
                for m = 1:2
                    b = model(m).modelrdv \ datardv;
                    % and obtain the best-fitting model RDM for each subject
                    fitted = model(m).modelrdv * b;
                    % evaluate CV prediction performance
                    for sub = 1:nsub
                        trainind = setdiff(1:nsub,sub);
                        % here we could average at one of two levels (estimates or
                        % predictions).
                        % For now, let's average predictions.
                        prediction = mean(fitted(:,trainind),2);
                        % dependent test case
                        rcv(sub,m) = atanh(corr(prediction,datardv(:,sub)));
                        % independent test case (subject never used in training
                        % split)
                        % (with fisher transform)
                        rvalid(sub,m) = atanh(corr(prediction,datardv(:,sub+nsub)));
                    end
                end
                % finally, model comparison by contrasting CV performances
                rcv(:,3) = rcv(:,1)-rcv(:,2);
                rvalid(:,3) = rvalid(:,1)-rvalid(:,2);
                sim_stdev_rcv(sigind1,sigind2,:) = std(rcv);
                sim_stdev_rvalid(sigind1,sigind2,:) = std(rvalid);
                sim_means_rcv(sigind1,sigind2,:) = mean(rcv);
                sim_means_rvalid(sigind1,sigind2,:) = mean(rvalid);
            end
        end
        % awkwardness courtesy of matlab parfor limitations
        means_rcv{simind} = sim_means_rcv;
        means_rvalid{simind} = sim_means_rvalid;
        stdev_rcv{simind} = sim_stdev_rcv;
        stdev_rvalid{simind} = sim_stdev_rvalid;
    end
    fprintf('finished in %.1f minutes (%.2fs per simulation)\n',toc(tstart)/60,...
        toc(tstart)/nsim);
    % now we unpack the cell arrays
    means_rcv = cat(4,means_rcv{:});
    means_rvalid = cat(4,means_rvalid{:});
    stdev_rcv = cat(4,stdev_rcv{:});
    stdev_rvalid = cat(4,stdev_rvalid{:});
    fprintf('saved results to cache at %s\n',cache);
    save(cache,'means_rcv','means_rvalid','stdev_rcv','stdev_rvalid');
else
    fprintf('using cached results\n');
    load(cache);
end

t_rcv = means_rcv ./ (stdev_rcv ./ sqrt(nsub));
t_rvalid = means_rvalid ./ (stdev_rvalid ./ sqrt(nsub));
p_rcv = tcdf(-t_rcv(:,:,1:2,:),nsub-1);
p_rvalid = tcdf(-t_rvalid(:,:,1:2,:),nsub-1);
% special handling of 2-tailed contrast case
p_rcv(:,:,3,:) = tcdf(-abs(t_rcv(:,:,3,:)),nsub-1)*2;
p_rvalid(:,:,3,:) = tcdf(-abs(t_rvalid(:,:,3,:)),nsub-1)*2;

prej_rcv = sum(p_rcv<.05,4) / nsim;
prej_rvalid = sum(p_rvalid<.05,4) / nsim;

% so for each parameter estimate, plot a colour map. For par one we should see a
% left-right gradient with effect size, with perhaps slightly biased estimates
% for the zero case on the left. For par two it will be up-down. And for
% par three (the contrast), we should expect unbiased estimates along the
% diagonal. Probably make this point by averaging the estimates for the zero
% signal case.

F = figure(100);
clf(F);
set(F,'defaultaxesfontsize',7,'defaulttextfontsize',7);
np = 0;
titles = {'model 1 effect','model 2 effect','contrast effect (1-2)'};
% first row
linind{1} = sub2ind([nsiglevel,nsiglevel],ones(1,nsiglevel),[1:nsiglevel]);
% first column
linind{2} = sub2ind([nsiglevel,nsiglevel],[1:nsiglevel],ones(1,nsiglevel));
% main diagonal
linind{3} = sub2ind([nsiglevel,nsiglevel],[1:nsiglevel],[1:nsiglevel]);
plotdims = [4,4];
% first plot legends
lims = [0 1];
difflim = .5 * [-1 1];
ncolor = 256;
np = np+1;
% this is a little loopy - we plot into the BOTTOM half of
% subplot(plotdims(1),plotdims(2),np) by pretending we have twice as many rows
% in the subplot grid
axl = subplot(plotdims(1)*2,plotdims(2),plotdims(1)+np);
imagesc(linspace(lims(1),lims(2),ncolor),0,linspace(lims(1),lims(2),ncolor));
set(axl,'tickdir','out','xtick',lims,'ytick',[],'plotboxaspectratio',[5 1 1]);
box(axl,'off');
title(axl,'rejection probability');
xlabel({'',''});
ylabel({'',''});

% then do the same on the third panel, and skip the fourth
np = np+2;
axl = subplot(plotdims(1)*2,plotdims(2),plotdims(1)+np);
imagesc(linspace(difflim(1),difflim(2),ncolor),0,...
    linspace(difflim(1),difflim(2),ncolor));
set(axl,'tickdir','out','xtick',[difflim(1) 0 difflim(2)],'ytick',[],...
    'plotboxaspectratio',[5 1 1]);
box(axl,'off');
title(axl,'difference');
xlabel({'',''});
ylabel({'',''});

% skip the fourth (over bars)
np = np+1;
for m = 1:3
    m_prej_rcv = prej_rcv(:,:,m);
    m_prej_rvalid = prej_rvalid(:,:,m);
    % column 1 - rejection probability for CV
    np = np+1;
    ax = subplot(plotdims(1),plotdims(2),np);
    h = plotpanel(ax,siglevel,m_prej_rcv,lims,linind{m});
    % if we on the first column, provide a row label
    ylabel({titles{m},'model 1 signal'});
    % if we are on the last row, provide column labels
    xlabel('');
    if np > (plotdims(1)-1)*plotdims(2)
        xlabel({'model 2 signal','dependent (CV)'});
    end

    % column 2 - validation estimate
    np = np+1;
    ax = subplot(plotdims(1),plotdims(2),np);
    h = plotpanel(ax,siglevel,m_prej_rvalid,lims,linind{m});
    % no y label in any case
    ylabel({'',''});
    % maybe x label if we are on the last row
    xlabel('');
    if np > (plotdims(1)-1)*plotdims(2)
        xlabel({'model 2 signal','independent (validation)'});
    end

    % column 3 - difference
    np = np+1;
    ax = subplot(plotdims(1),plotdims(2),np);
    h = plotpanel(ax,siglevel,m_prej_rcv-m_prej_rvalid,difflim,...
        linind{m});
    % no y label in any case
    ylabel({'',''});
    % maybe x label if we are on the last row
    xlabel('');
    if np > (plotdims(1)-1)*plotdims(2)
        xlabel({'model 2 signal','bias (dependent-independent)'});
    end

    % column 4 - summary stats of null
    np = np+1;
    ax = subplot(plotdims(1),plotdims(2),np);
    alpha(1) = mean(m_prej_rcv(linind{m}));
    alpha(2) = mean(m_prej_rvalid(linind{m}));
    hold(ax,'on');
    bh(1) = bar(ax,1,alpha(1),'facecolor',[0 .5 0],'edgecolor','none');
    bh(2) = bar(ax,2,alpha(2),'facecolor',[0 0 .5],'edgecolor','none');
    ylabel('false positive rate');
    set(ax,'tickdir','out','xtick',[1 2],'xticklabel',{'CV','valid'},...
        'plotboxaspectratio',[1 1 1]);
    box(ax,'off');
    xlim(ax,[.5 2.5]);
    lh = line(xlim,[.05 .05],'color','k','linestyle',':');
    th = text(max(xlim) + range(xlim)*.1,.05,'alpha');
    xlabel(ax,{'',''});
    if np == plotdims(1)*plotdims(2)
        xlabel({'means under H0',''});
    end
end
colormap(jet(ncolor));
% make sure print -dpdf will work ok
figsize = 15 * [1.2 1];
set(F,'paperunits','centimeters','paperposition',[0 0 figsize],...
    'papersize',figsize,'name','test_cvdependence: full result');
print(fullfile(fundir,'test_cvdependence_full.pdf'),'-dpdf','-r300');

F2 = figure(200);
figsize2 = figsize(2) * [.4 .4];
clf(F2);
set(F2,'defaultaxesfontsize',7,'defaulttextfontsize',7);
ax = gca(F2);
plot(prej_rcv(:),prej_rvalid(:),'.','markersize',8);
axis([0 1 0 1]);
lh = line([0 1],[0 1],'color','k','linestyle',':');
uistack(lh,'bottom');
box(ax,'off');
set(ax,'tickdir','out','plotboxaspectratio',[1 1 1],...
    'dataaspectratio',[1 1 1],'xtick',[0 .25 .5 .75 1],'ytick',[0 .25 .5 .75 1]);
xlabel(ax,{'rejection probability','dependendent (cv)'});
ylabel(ax,{'rejection probability','independent (validation)'});
set(F2,'paperunits','centimeters','paperposition',[0 0 figsize2],...
    'papersize',figsize2,...
    'name','test_cvdependence: rejection probabilities');
print(fullfile(fundir,'test_cvdependence_rejprob.pdf'),'-dpdf','-r300');

if debug
    keyboard;
end

% support function for the imagesc plots
function h = plotpanel(ax,siglevel,y,lims,linind)
h = imagesc(siglevel,siglevel,y,lims);
hold(ax,'on');
[x,y] = meshgrid(siglevel(:),siglevel(:));
p = plot(x(linind),y(linind),'ow');
set(ax,'dataaspectratio',[1 1 1],'plotboxaspectratio',[1 1 1],...
    'tickdir','out','ydir','normal');
box(ax,'off');
