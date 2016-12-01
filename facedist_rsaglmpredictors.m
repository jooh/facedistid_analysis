% predictors for multirsa on facedistid data.
%
% rdmgroups = facedist_rsaglmpredictors(rdms,ismeanrsa)
function rdmgroups = facedist_rsaglmpredictors(rdms,ismeanrsa)

if ieNotDefined('ismeanrsa')
    ismeanrsa = false;
end

ncon = size(rdms(1).RDM,1);
switch ncon
    case 24
        fullrdm = true;
    case 12
        fullrdm = false;
    otherwise
        error('unknown rdm size: %d',size(rdms(1).RDM,1));
end

eccrdm = rdms(strcmp({rdms.name},'no_direction_tuning'));
eccrdm.name = 'eccentricity';
angrdm = rdms(strcmp({rdms.name},'no_eccentricity_discrimination'));
angrdm.name = 'direction';

% ensure the sub-typical and typical-caricature pairs are distinct
indsa = 1:4;
indsb = 5:8;
% smallest floating point number that will still count as different
f = eps(max(asrdmvec(eccrdm)));
eccrdm.RDM(indsb,indsa) = eccrdm.RDM(indsb,indsa)+f;
eccrdm.RDM(indsa,indsb) = eccrdm.RDM(indsa,indsb)+f;
if fullrdm
    % need to do the remaining quadrants too
    eccrdm.RDM(indsb+12,indsa+12) = eccrdm.RDM(indsb+12,indsa+12)+f;
    eccrdm.RDM(indsa+12,indsb+12) = eccrdm.RDM(indsa+12,indsb+12)+f;
    eccrdm.RDM(indsb,indsa+12) = eccrdm.RDM(indsb,indsa+12)+f;
    eccrdm.RDM(indsa+12,indsb) = eccrdm.RDM(indsa+12,indsb)+f;
    eccrdm.RDM(indsb+12,indsa) = eccrdm.RDM(indsb+12,indsa)+f;
    eccrdm.RDM(indsa,indsb+12) = eccrdm.RDM(indsa,indsb+12)+f;
end

eccfir = rdm2fir(eccrdm,[],true);
% resort
eccfir = eccfir([1 2 6 4 3 5]);
% rename
newnames = {'sub-sub','typical-typical','caricature-caricature',...
    'sub-typical','typical-caricature','sub-caricature'};
[eccfir.name] = deal(newnames{:});

angfir = rdm2fir(angrdm,'%02.0f',true);
rdmgroups = emptystruct('name','RDM');

% make view FIR RDMs
if fullrdm
    viewrdm(1) = rdms(strcmp({rdms.name},'view'));
    viewrdm(1).name = 'view_across';
    viewrdm(2) = struct('name','view_left','RDM',...
        [ones(12,12),zeros(12,12); zeros(12,24)]);
    viewrdm(2).RDM(diagind(size(viewrdm(2).RDM))) = 0;
    viewrdm(3) = struct('name','view_right','RDM',...
        [zeros(12,24); zeros(12,12),ones(12,12)]);
    viewrdm(3).RDM(diagind(size(viewrdm(3).RDM))) = 0;
    rdmgroups(end+1).name = 'view';
    rdmgroups(end).RDM = viewrdm;
end

rdmgroups(end+1).name = 'eccentricity';
rdmgroups(end).RDM = eccfir;
rdmgroups(end+1).name = 'direction';
rdmgroups(end).RDM = angfir;

% combinations of 2 predictors
if fullrdm
    rdmgroups(end+1).name = 'view*eccentricity';
    rdmgroups(end).RDM = rdmsplitbyfir(eccfir,viewrdm);
    rdmgroups(end+1).name = 'view*direction';
    rdmgroups(end).RDM = rdmsplitbyfir(angfir,viewrdm);
end
rdmgroups(end+1).name = 'eccentricity*direction';
rdmgroups(end).RDM = rdmsplitbyfir(angfir,eccfir);

% and having achieved this, I believe we can now go all the way and build
% the full model
if fullrdm
    rdmgroups(end+1).name = 'view*eccentricity*direction';
    rdmgroups(end).RDM = rdmsplitbyfir(angfir,rdmgroups(...
        strcmp({rdmgroups.name},'view*eccentricity')).RDM);
end

if ismeanrsa && fullrdm
    offdiagind = circshift(diagind([24 24]),[12 0]);
    for g = 1:numel(rdmgroups)
        if any(strfind(rdmgroups(g).name,'direction'))
            % leave these alone since we do want to estimate 0 direction
            % case
            continue
        end
        for r = 1:numel(rdmgroups(g).RDM)
            rdmgroups(g).RDM(r).RDM(offdiagind) = 0;
        end
    end
end

% also do a full, unaveraged model (this is kind of like the noise ceiling
% or 'no model' case). So when you use cvpredictruns with this model you
% effectively just correlate the data RDMs between the splits.
npairs = nchoosek(ncon,2);
base = eye(npairs);
% convert to rdm (bit hacky here to prevent this from being interpreted as
% a single ncon by ncon RDM rather than ncon vectorised RDMs.
baserdm = arrayfun(@(x)asrdmmat(base(:,x)),1:npairs,'uniformoutput',0);
rdmgroups(end+1).name = 'reproducibility';
rdmgroups(end).RDM = struct('name',mat2strcell(1:npairs,...
    'reproducibility_%03d'),'RDM',baserdm);

% finally sure everything is single precision for speed, and reinstate mean
% dissimilarity data (likely got lost in the process)
fullind = strcmp({rdms.name},'full');
for g = 1:numel(rdmgroups)
    for r = 1:numel(rdmgroups(g).RDM)
        rdmgroups(g).RDM(r).RDM = single(rdmgroups(g).RDM(r).RDM);
        inds = rdmgroups(g).RDM(r).RDM ~= 0;
        m = mean(rdms(fullind).RDM(inds));
        % 0 predictions are a little hard to handle. We fudge this by
        % setting it to a small but non-zero value
        if m==0
            m = f;
        end
        rdmgroups(g).RDM(r).RDM(inds) = m;
        assert(~isnan(m));
        % also add the mean dissimilarity from the euclidean RDM to the
        % file name. This is a bit hacky (really the model RDM should have
        % a covariate field or something like that) but for now...
        rdmgroups(g).RDM(r).name = sprintf('%s_m%.1f',...
            rdmgroups(g).RDM(r).name,m);
    end
end

if fullrdm && ismeanrsa
    for g = 1:numel(rdmgroups)
        % make sure off entries are NaN rather than 0
        for r = 1:numel(rdmgroups(g).RDM)
            rdmgroups(g).RDM(r).RDM(rdmgroups(g).RDM(r).RDM==0) = NaN;
            rdmgroups(g).RDM(r).RDM = zerodiagonal(rdmgroups(g).RDM(r).RDM);
        end
        % and collapse left / right into a 'within' predictor
        leftind = strfindcell({rdmgroups(g).RDM.name},'left');
        rightind = [];
        % NB this loop only enters if there are leftinds to process
        for l = 1:numel(leftind)
            left = rdmgroups(g).RDM(leftind(l));
            [~,rightind(l)] = intersect({rdmgroups(g).RDM.name},strrep(left.name,...
                'left','right'));
            right = rdmgroups(g).RDM(rightind(l));
            newrdm = NaN([24 24],class(left.RDM));
            % to test
            newrdm(1:12,1:12) = left.RDM(1:12,1:12);
            newrdm(13:24,13:24) = right.RDM(13:24,13:24);
            rdmgroups(g).RDM(end+1) = struct('name',strrep(left.name,'left',...
                'within'),'RDM',newrdm);
        end
        % get rid of left and right to keep the number of comparisons below
        % the triple digits (we can remove this step here if we ever do want to
        % look at unilateral effects. But let's not).
        rdmgroups(g).RDM([leftind;rightind]) = [];
    end
end

% add image of each for ease of visualisation later
allcolors = cmap_glasbey;
for g = 1:numel(rdmgroups)
    if strcmp(rdmgroups(g).name,'reproducibility')
        continue;
    end
    thisrdmat = single(asrdmmat(rdmgroups(g).RDM)~=0);
    nmat = size(thisrdmat,3);
    thisc = allcolors(1:(nmat+1),:);
    allcolors(1:(nmat+1),:) = [];
    % add offset for each
    urdmat = bsxfun(@times,thisrdmat,reshape(1:nmat,[1 1 nmat]));
    % and sum to obtain final matrix
    fullmat = sum(urdmat,3);
    % ensure diagonal is 0
    thisc(1,:) = 1;
    % dump as image
    rdmgroups(g).im = asrdmimage(fullmat,thisc);
end
