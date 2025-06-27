function paretoFront = GetParetoFront(non_dom_pop)
    N = size(non_dom_pop,1);
    isDominated = false(1, N);

    for i = 1:N
        for j = 1:N
            if i == j, continue; end

            if Dominates(non_dom_pop(j), non_dom_pop(i)) == 1
                % if j dominate i
                isDominated(i) = true;
                break;  
            end
        end
    end

    paretoFrontAll = non_dom_pop(~isDominated);
    
    %% delete pop that have same Cost
    costs = reshape(vertcat(paretoFrontAll.Cost),[2 size(paretoFrontAll,1)])';
    [~, uniqueIdx] = unique(costs, 'rows', 'stable');
    paretoFront = paretoFrontAll(uniqueIdx);
end
