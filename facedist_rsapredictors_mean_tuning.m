% predictors for MeanRSA on facedistid data. Uses facedist_rsaglmpredictors
% internally to generate the correct FIR-like model RDMs.
%
% rdmgroups = facedist_rsapredictors_mean_tuning(rdms)
function outg = facedist_rsapredictors_mean_tuning(rdms)

switch size(rdms(1).RDM,1)
    case 24
        fullrdm = true;
    case 12
        fullrdm = false;
    otherwise
        error('unknown rdm size: %d',size(rdms(1).RDM,1));
end

rdmgroups = facedist_rsaglmpredictors(rdms,true);

% combine the relevant groups into one big set of
% predictors for meanrsa
if fullrdm
    % so we want everything view specific, so in all cases we want view*
    % interactions. So it's going to be view*eccentricity and
    % view*eccentricity*direction. Too many RDMs perhaps (111). In which case we
    % drop left and right.
    % but first set across view diagonal to NaN so we can do a like for
    % like comparison within and across views.
    ecc = rdmgroups(4).RDM;
    % indices to off-diagonals
    %offdiag = circshift(diagind([24 24]),[12 0]);
    %for e = 1:numel(ecc)
        %ecc(e).RDM(offdiag) = NaN;
    %end
    outg = [ecc rdmgroups(7).RDM];
else
    % take the eccentricity blocks and the direction*eccentricity
    % interaction (no global direction RDM since this has a problematic 0
    % case (only across eccentricity levels).
    outg = [rdmgroups([1 3]).RDM];
end
