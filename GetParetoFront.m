function paretoFront = GetParetoFront(non_dom_pop)
    N = size(non_dom_pop,1);
    isDominated = false(1, N);

    for i = 1:N
        for j = 1:N
            if i == j, continue; end

            if checkDomination(non_dom_pop(j).Cost, non_dom_pop(i).Cost) == 1
                % if j dominate i
                isDominated(i) = true;
                break;  
            end
        end
    end

    paretoFrontAll = non_dom_pop(~isDominated);
    
    % delete pop that have same Cost
    costs = vertcat(paretoFrontAll.Cost);
    [~, uniqueIdx] = unique(costs, 'rows', 'stable');
    paretoFront = paretoFrontAll(uniqueIdx);
end
