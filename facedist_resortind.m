% indices for rearranging facedist RDMs according to view (if present) and
% then eccentricity. n must be 12 or 24.
%
% resortind = facedist_resortind(n)
function resortind = facedist_resortind(n)

switch n
    case 12
        % sort according to eccentricity
        [~,resortind] = sort(repmat(1:3,[1 4]));
    case 24
        % sort according to view, then eccentricity
        [~,resortind] = sort([repmat(1:3,[1 4]) repmat(4:6,[1 4])]);
    otherwise
        error('unsupported n: %d',n);
end
