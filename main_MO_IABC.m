%% MO_IABC traditional MOO using domination sort 
%% DEPLOYMENT

%%
clc;
clear;

%% 
close all;
%% Network parameter
%for N=20:5:30
%for trial=1:10

% Monitor area
Covered_Area = zeros(100,100);
%Obstacle_Area = gen_random_distribution_area(30,29);
Obstacle_Area = ones(100,100);
%Obstacle_Area=genarea();
%Obstacle_Area = gen_target_area(1000);

% nodes info
MaxIt = 100;              % Maximum Number of Iterations
a = 1;                    % Acceleration Coefficient Upper Bound
N = 60;
nPop=50;
rc = 10;
sink=[50 50];
%name=['./data/hetero_target/',num2str(N),' nodes/case4/hetero_',num2str(trial),'.mat'];
%name=['./data/comparison/barrier/ABC/ABC',num2str(N),'nodes_',num2str(trial),'.mat'];

% homogenerous sensor 
rs=ones(1,N)*10;

% bee parameter
Scout_bee=nPop;             % number of scout bees
Onlooker_bee=nPop;          % number of onlooker bees
BestCost=zeros(MaxIt,1);    % Best solution every iteration
L=zeros(1,nPop);            % depletion matrix

%% Init first pop
empty_individual.Position=[]; % Empty solution
empty_individual.Cost=[];       % Empty cost function of that solution
pop=repmat(empty_individual,nPop,1); % pop includes nPop solutions
non_dom_pop=repmat(empty_individual,1,1); % non_dom_pop includes non-dominated solutions
clear empty_individual;

for k = 1:nPop
    alpop = unifrnd(sink(1)-rc/2,sink(2)+rc/2,[N 2]);
    % gen sink node
    alpop (1,:)= sink;
    pop(k).Position=alpop;
    pop(k).Cost(1)=Cov_Func_v2(pop(k).Position,rs,Obstacle_Area,Covered_Area);
    pop(k).Cost(2) = Energy_consumption(pop(k).Position,rc);
    pop(k).Cost(3) = Life_Time(pop(k).Position,rc);
end
alpop = unifrnd(sink(1)-rc/2,sink(2)+rc/2,[N 2]);
alpop (1,:)= sink;
non_dom_pop(1).Position=alpop;
non_dom_pop(1).Cost(1)=Cov_Func_v2(non_dom_pop(1).Position,rs,Obstacle_Area,Covered_Area);
non_dom_pop(1).Cost(2) = Energy_consumption(non_dom_pop(1).Position,rc);
non_dom_pop(1).Cost(3) = Life_Time(non_dom_pop(1).Position,rc);
clear k alpop sink;
%%     ABC Main Loop
for it = 1:MaxIt
    %% Global search
    for i=1:Scout_bee
        k=randi([1,nPop]);
        phi=a*unifrnd(-1, +1, [N 2])*(1-L(i)/MaxIt)^5;
        alpop = pop(i).Position + phi.*(pop(i).Position-pop(k).Position);
        alpop(:,1:2) = min(max(alpop(:,1:2), 1),size(Obstacle_Area,1)-1);
        if Connectivity_graph(Graph(alpop(:,1:2),rc),[]) == 1
            alpop_Cost(1) = Cov_Func_v2(alpop,rs,Obstacle_Area,Covered_Area);
            alpop_Cost(2) = Energy_consumption(alpop,rc);
            alpop_Cost(3) = Life_Time(alpop,rc);

            if checkDomination(alpop_Cost,pop(i).Cost) == 1
                % alpop dominate pop
                pop(i).Position=alpop;
                pop(i).Cost=alpop_Cost;
            elseif checkDomination(alpop_Cost,pop(i).Cost) == -1
                % pop dominate alpop
                L(i)=L(i)+1;
                continue;
            else
                % alpop nondominate pop
                non_dom_pop(end+1,1).Position=alpop;
                non_dom_pop(end,1).Cost=alpop_Cost;
                non_dom_pop = GetParetoFront(non_dom_pop);        
            end
        else
            L(i)=L(i)+1;
        end
    end

    clear alpop_cov alpop alpop2 k phi;
    %% ranking pop
    E=zeros(1,nPop);
    for i=1:nPop
        E(i)=sum(pop(i).Cost);
    end
    E_ave=E/sum(E);

    %% Local search
    for j=1:Onlooker_bee
        % randomly choose a pop, prioritize that have high coverage
        if all(isnan(E_ave))
            i=randi([1,nPop]);
        else
            i=randsrc(1,1,[1:nPop;E_ave]);
        end
        % change the position of each solution
        for k=1:N
            alpop=pop(i).Position;
            h=randi([1,N]);
            phi=a*unifrnd(-1, +1, [1 2])*(1-L(i)/MaxIt)^2;
            alpop(k,1:2)  = pop(i).Position(k,1:2) + phi.*(pop(i).Position(k,1:2)-pop(i).Position(h,1:2));
            alpop(:,1:2) = min(max(alpop(:,1:2), 1),size(Obstacle_Area,1)-1);
            if Connectivity_graph(Graph(alpop(:,1:2),rc),[]) == 1
                alpop_Cost(1)=Cov_Func_v2(alpop,rs,Obstacle_Area,Covered_Area);
                alpop_Cost(2) = Energy_consumption(alpop,rc);
                alpop_Cost(3) = Life_Time(alpop,rc);
    
                if checkDomination(alpop_Cost,pop(i).Cost) == 1
                    % alpop dominate pop
                    pop(i).Position=alpop;
                    pop(i).Cost=alpop_Cost;
                elseif checkDomination(alpop_Cost,pop(i).Cost) == -1
                    % pop dominate alpop
                    continue;
                else
                    % alpop nondominate pop
                    non_dom_pop(end+1,1).Position=alpop;
                    non_dom_pop(end,1).Cost=alpop_Cost;
                    non_dom_pop = GetParetoFront(non_dom_pop);
                end
            end
        end
    end
    clear alpop_cov alpop alpop2 k h phi check_dom h j k phi E E_ave;
    disp(['Iteration ' num2str(it) ':' num2str(size(non_dom_pop,1)) ' Non-dominated solutions']);  
    
    %% plot pareto
    data=zeros(size(non_dom_pop,1),3);
    for i = 1: size(non_dom_pop,1)
        data(i,1) = non_dom_pop(i).Cost(1);
        data(i,2) = non_dom_pop(i).Cost(2);
        data(i,3) = non_dom_pop(i).Cost(3);
    end
    plot3(data(:,1),data(:,2),data(:,3)*100,'o','Color','r')
    %plot(data(:,1),data(:,2),'o','Color','r');
    xlabel('Non-coverage cost function') 
    ylabel('Energy cost function') 
    zlabel('Life time')
    drawnow;
end
%% plot
% % %alpop=pop(4).Position;
% alpop=BestSol.Position;
% [coverage,Covered_Area] = Cov_Func(alpop,rs,theta0,Obstacle_Area,Covered_Area);
% % show map
% %imagesc(imresize(Obstacle_Area,[100 100],'bilinear'))
% % figure;
% % plot(BestCost);
% figure;
% hold on;
% 
% % show map
% 
% %[obs_row, obs_col] = find(Covered_Area == 0);
% %plot(obs_row-1, obs_col-1,'.', 'MarkerSize', 2, 'Color', 'green');
% 
% 
% %[discovered_obs_row, discovered_obs_col] = find(Covered_Area == -1);                    % show discovered map
% %plot(discovered_obs_row-1, discovered_obs_col-1,'.', 'MarkerSize', 20, 'Color', 'red');
% %colorbar;
% %{
% for i = 1:1:numel(G.Edges.EndNodes)/2
%     plot([pop(G.Edges.EndNodes(i,1)*2-1),pop(G.Edges.EndNodes(i,2)*2-1)],[pop(G.Edges.EndNodes(i,1)*2),pop(G.Edges.EndNodes(i,2)*2)],'Color','blue','linewidth',1);
% end
% %}
% [obs_row, obs_col] = find(Obstacle_Area == 1);
% plot(obs_row, obs_col,'.', 'MarkerSize', 5, 'Color', 'red');
% [obs_row, obs_col] = find(Covered_Area == 1);
% plot(obs_row, obs_col,'.', 'MarkerSize', 5, 'Color', 'blue');
% clear i x_fill y_fill theta obs_col obs_row;
% 
% axis equal;
% xlim([0 size(Obstacle_Area,1)]);
% ylim([0 size(Obstacle_Area,2)]);
% title(['Coverage ratio: ' num2str(BestSol.Cost*100) '% ' ]);
% grid on;
% drawnow;
% 
% % save fig
% exportgraphics(gcf, './data/zFigures/deployment figures/hetero_target_case4.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');

%%
%save(name)
%end
%end
figure;
%hold on;


