%% Directional Sensor
%% DEPLOYMENT

%%
clc;
clear;

%% 
close all;
%% Network parameter
% for trial=2:10
% name=['./data/comparison/GA/GA30nodes_',num2str(trial),'.mat' ];
% Monitor area
Covered_Area = zeros(100,100);
%Obstacle_Area = gen_random_distribution_area(30,29);
Obstacle_Area = ones(100,100);

% nodes info
MaxIt = 1000;              % Maximum Number of Iterations
N = 60;
nPop = 200;
rc = 30; 
rs = ones(1,N)*10;
stat(1,:)=rs;
stat(2,1)=rc;
%name=['./data/hetero_target/',num2str(N),' nodes/case4/hetero_',num2str(trial),'.mat'];

sink=[50 50];

% GA parameter
crossover_rate=0.9;
mutation_rate=0.1;

%% Init first pop
empty_individual.Position=[]; % Empty solution
empty_individual.Cost=[];     % Empty cost function of that solution
pop=repmat(empty_individual,nPop,1); % pop includes nPop solutions
non_dom_pop=repmat(empty_individual,nPop,1); % non_dom_pop includes non-dominated solutions
clear empty_individual;

for k = 1:nPop
    alpop = unifrnd(20,80,[N 2]);
    % gen sink node
    alpop (1,:)= sink;
    pop(k).Position = alpop;
    G = Graph(alpop,rc);
    pop(k).Cost(1) = Cov_Func_v2(pop(k).Position,rs,Obstacle_Area,Covered_Area);
    pop(k).Cost(2) = Energy_consumption(G);
    pop(k).Cost(3) = Life_Time(G);
    non_dom_pop(k).Position = alpop;
    non_dom_pop(k).Cost = pop(k).Cost;
end
%%
non_dom_pop = GetParetoFront(non_dom_pop);
clear k alpop sink;
%%     GA Main Loop
for it = 1:MaxIt
    for i=1:nPop
        alpop=pop(i).Position;
        for j=1:N
            if rand()<crossover_rate
                l=randi([1,nPop]);
                alpop(j,:) = pop(l).Position(j,:) ;
            end
            if rand()<mutation_rate
                l=randi([1,nPop]);
                h=randi([1,N]);
                alpop(j,:) = pop(l).Position(h,:) ;
            end
        end
        alpop(:,1:2) = min(max(alpop(:,1:2), 1),size(Obstacle_Area,1)-1);
        if Connectivity_graph(Graph(alpop(:,1:2),rc),[]) == 1
            G = Graph(alpop,rc);
            alpop_Cost(1) = Cov_Func_v2(alpop,rs,Obstacle_Area,Covered_Area);
            alpop_Cost(2) = Energy_consumption(G);
            alpop_Cost(3) = Life_Time(G);
            
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

    clear alpop_cov alpop alpop2 k phi check_dom;
    disp(['Iteration ' num2str(it) ':' num2str(size(non_dom_pop,1)) ' Non-dominated solutions']);
    clear h j k phi alpop alpop2 alpop_cov E E_ave;
    %%
    %figure;
    %hold on;
    data=zeros(size(non_dom_pop,1),3);
    for i = 1: size(non_dom_pop,1)
        data(i,1) = non_dom_pop(i).Cost(1);
        data(i,2) = non_dom_pop(i).Cost(2);
        data(i,3) = non_dom_pop(i).Cost(3);
    end
    plot3(data(:,1),data(:,2),data(:,3),'o','Color','r')
    %plot(data(:,1),data(:,2),'o','Color','r');
    xlabel('Non-coverage cost function') 
    ylabel('Energy cost function') 
    zlabel('Life time cost function')
    drawnow;
end
%% plot
% %alpop=pop(4).Position;
% alpop=BestSol.Position;
% coverage = Cost_Func(alpop,weight,stat,Obstacle_Area,Covered_Area);
% % show map
% %imagesc(imresize(Obstacle_Area,[100 100],'bilinear'))
% figure;
% plot(BestCost);
% figure;
% hold on;
% 
% % show map
% [obs_row, obs_col] = find(Obstacle_Area == 1);
% plot(obs_row, obs_col,'.', 'MarkerSize', 10, 'Color', 'blue');
% [obs_row, obs_col] = find(Covered_Area == 1/2);
% plot(obs_row, obs_col,'.', 'MarkerSize', 2, 'Color', 'green');
% [obs_row, obs_col] = find(Covered_Area == 1);
% plot(obs_row, obs_col,'.', 'MarkerSize', 10, 'Color', 'cyan');
% %[discovered_obs_row, discovered_obs_col] = find(Covered_Area == -1);                    % show discovered map
% %plot(discovered_obs_row-1, discovered_obs_col-1,'.', 'MarkerSize', 20, 'Color', 'red');
% %colorbar;
% %{
% for i = 1:1:numel(G.Edges.EndNodes)/2
%     plot([pop(G.Edges.EndNodes(i,1)*2-1),pop(G.Edges.EndNodes(i,2)*2-1)],[pop(G.Edges.EndNodes(i,1)*2),pop(G.Edges.EndNodes(i,2)*2)],'Color','blue','linewidth',1);
% end
% %}
% 
% 
% %clear i x_fill y_fill theta obs_col obs_row;
% 
% axis equal;
% xlim([0 size(Obstacle_Area,1)]);
% ylim([0 size(Obstacle_Area,2)]);
% title(['Weighted Coverage ratio:' num2str(BestSol.Cost*100) '% ' ]);
% grid on;
% drawnow;
% %{
% if save_fig==1
%     saveas(gcf, ['Iteration' num2str(it) '.pdf' ]); % Lưu dưới dạng PDF
% end
% %}
%%
%save(name)
%end

