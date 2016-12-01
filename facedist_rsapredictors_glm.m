% predictors for GLM RSA on facedistid data. 
%
% outg = facedist_rsapredictors_glm(rdms)
function outg = facedist_rsapredictors_glm(rdms)

switch size(rdms(1).RDM,1)
    case 24
        fullrdm = true;
    case 12
        fullrdm = false;
    otherwise
        error('unknown rdm size: %d',size(rdms(1).RDM,1));
end

% get the components
target = rdms(strcmp({rdms.name},'full')).RDM.^2;
% target with nan diagonal so we can average non-zero values conveniently
% with nanmean
nantarget = target;
nantarget(diagind(size(nantarget))) = NaN;

meanoff = @(x)mean(x(logical(1-eye(4))));

% delta ecc - min discriminability at each ecc level (effectively the 0
% direction change case since this is present in each block)
deltaecc = blkfun(target,[4 4],@min);

% get whatever variance is left in the full face space (these are all
% direction effects basically).
abseccdeltadir = target - deltaecc;

constant = zerodiagonal(ones(size(target)));

splitdireffects = false;
if splitdireffects
    % Now to go a step further, we think that there may be  an increase in
    % discriminability with eccentricity, but not necessarily a direction
    % _tuning_ effect, because we may be outside the saturation point of this
    % function.

    % so to capture such an effect, we simply block average out the direction
    % effect 
    % eccabs = zerodiagonal(blkfun(abseccdeltadir,[4 4],@mean));
    % let's try
    eccabs = zerodiagonal(blkfun(abseccdeltadir,[4 4],meanoff));
    % (indeed, these two are identical apart from a scale factor)
    % and to complete the variance partition, subtract off this mean
    % eccentricity effect to get the variance that's uniquely attributable to
    % delta(dir) without effects of absolute or delta ecc.
    deltadir = zerodiagonal(abseccdeltadir - eccabs);
    % All this raises the question of whether we should include the main
    % effect of direction as well. Unfortunately this won't partition
    % neatly under the above scheme. Anyway, it doesn't make sense since
    % delta dir is a polar predictor. It has to interact with radius in a
    % cartesian space - imagine the case where absecc=0. How could there be
    % an effect of deltadir? This is why the coef angdist model was always
    % a bit silly.

    % now deltadir seems a bit hacky in that we get there through 2 previous
    % steps. We can simplify by just subtracting the mean for each block, which
    % turns out to give an identical answer (down to rounding error)
    deltadir2 = target - zerodiagonal(blkfun(target,[4 4],meanoff));

    deltadirmain = repmat(asrdmmat((asrdmvec(target(1:4,1:4)))),...
        [3*[fullrdm+1] 3*[fullrdm+1]]);
    % so this is effectively the dir x absecc interaction term
    abseccdeltadir2 =  deltadirmain .* eccabs;
    % And a more complete parameterisation with both main effects and
    % absecc*delta dir interactions would be
    predmatnew = sqrtsigned(cat(3,deltaecc,eccabs,deltadirmain,...
        abseccdeltadir2,constant));
    % Both models can perfectly reconstruct the original model distances.
    % The problem with this parameterisation is that the interaction term ends
    % up accounting for all the absecc and deltadir effects. This is a
    % necessary consequence of course since one is basically direction and the
    % other radius.

    % When we fit facepairs with each we get weird negative estimates for delta
    % dir main with the new parameterisation. The original parameterisation
    % looks much better - about 50% over-representation of delta ecc, and the
    % other two making a similar contribution. So I think for interpretability
    % this is the one to go with.
    % and the 'deltadir' predictor isn't quite this (because absecc has been
    % subtracted off).
    % so we can imagine 2 basic models. The old model goes like this.
    predmat = sqrtsigned(cat(3,deltaecc,eccabs,deltadir,constant));
    names = {'delta_eccentricity','abs_eccentricity','direction','constant'};
else
    predmat = sqrtsigned(cat(3,deltaecc,abseccdeltadir,constant));
    names = {'delta_eccentricity','delta_direction','constant'};
end

% but here it no longer is...
% concatenate, return to data units

if fullrdm
    % mask out each view in turn
    % Actually, let's regularise by collapsing left and right
    f = false([12,12]);
    t = true([12,12]);
    withinmask = repmat([f t; t f],[1 1 size(predmat,3)]);
    %leftmask = true([24 24]);
    %leftmask(1:12,1:12) = false;
    %leftmask = repmat(leftmask,[1 1 size(predmat,3)]);
    %predmat_left = predmat;
    %predmat_left(leftmask) = 0;
    %rightmask = circshift(leftmask,[12 12]);
    %predmat_right = predmat;
    %predmat_right(rightmask) = 0;

    predmat_within = predmat;
    predmat_within(withinmask) = 0;
    predmat_across = predmat;
    predmat_across(~withinmask) = 0;

    % now the final design matrix is
    predmat= cat(3,predmat_within,predmat_across);
    prefit = @(x,p)cellfun(@(thisx)[p thisx],x,'uniformoutput',0);
    % and the names are
    names = [prefit(names,'view_within_') ...
        prefit(names,'view_across_')];
end

% convert to struct and we're done
outg = rdm2struct(predmat,names);
