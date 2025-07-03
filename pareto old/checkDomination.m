function result = checkDomination(f1, f2)
% CHECK DOMINATION
% result:
%   = 1  if f1 dominate f2
%   = -1 if f2 dominate f1
%   = 0  if non-dominated


    if all(f2==f1)    % f1 = f2
        result = 2;
    elseif all(f1 <= f2) && any(f1 < f2)
        result = 1;   % f1 dominates f2
    elseif (all(f2 <= f1) && any(f2 < f1))
        result = -1;  % f2 dominates f1
    else
        result = 0;   % non-dominated
    end
end
