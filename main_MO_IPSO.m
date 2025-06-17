%% Directional Sensor
%% DEPLOYMENT

%%
clc;
clear;

%% 
close all;
%% Network parameter
%for trial=1:10
%name=['./data/comparison/PSO/PSO25nodes_',num2str(trial),'.mat' ];
% Monitor area
Covered_Area = zeros(100,100);
%Obstacle_Area = gen_random_distribution_area(30,29);
Obstacle_Area = ones(100,100);

% nodes info
MaxIt = 100;              % Maximum Number of Iterations
N = 60;
nPop=50;
rc = 10;
rs=ones(1,N)*10;
stat(1,:)=rs;
stat(2,1)=rc;
weight=[0.5 0.5];

sink=[50 50];

% FOA parameter
w = 0.7;                  % Trọng số quán tính
wdamp=0.7;
c1 = 1;                 % Hệ số học hỏi cá nhân
c2 = 0.5;                 % Hệ số học hỏi xã hội

%% Init first pop
empty_individual.Position=[]; % Empty solution
empty_individual.Cost=[];     % Empty cost function of that solution
empty_individual.Velocity=[]; % Empty svelocity
pop=repmat(empty_individual,nPop,1); % pop includes nPop solutions
non_dom_pop=repmat(empty_individual,1,1); % non_dom_pop includes non-dominated solutions

clear empty_individual;

for k = 1:nPop
    alpop = unifrnd(sink(1)-rc,sink(2)+rc,[N 2]);
    % gen sink node
    alpop (1,:)= sink;
    pop(k).Position=alpop;
    pop(k).Cost(1)=Cov_Func_v2(pop(k).Position,rs,Obstacle_Area,Covered_Area);
    pop(k).Cost(2) = Energy_consumption(pop(k).Position,rc);
    pop(k).Cost(3) = Late_Func(pop(k).Position,rc);
    pop(k).Velocity=zeros(N,3);
end
alpop = unifrnd(sink(1)-rc/2,sink(2)+rc/2,[N 2]);
alpop (1,:)= sink;
non_dom_pop(1).Position=alpop;
non_dom_pop(1).Cost(1)=Cov_Func_v2(non_dom_pop(1).Position,rs,Obstacle_Area,Covered_Area);
non_dom_pop(1).Cost(2) = Energy_consumption(non_dom_pop(1).Position,rc);
non_dom_pop(1).Cost(3) = Late_Func(non_dom_pop(1).Position,rc);
clear k alpop sink;

%%     GA Main Loop
for it = 1:MaxIt
    for i=1:nPop
        for j=1:N
            alpop=pop(i).Position;
            k=randi([1,N]);
            h=randi([1,size(non_dom_pop,1)]);
            alpop_v=w*pop(i).Velocity(j,1:2) + c1*rand()*(pop(i).Position(k,1:2) - pop(i).Position(j,1:2)) + c2*rand()*(non_dom_pop(h).Position(j,1:2) - pop(i).Position(j,1:2));
            alpop(j,1:2) = pop(i).Position(j,1:2) + alpop_v;
            alpop(:,1:2) = min(max(alpop(:,1:2), 1),size(Obstacle_Area,1)-1);
            if Connectivity_graph(Graph(alpop(:,1:2),rc),[]) == 1
                alpop_Cost(1) = Cov_Func_v2(alpop,rs,Obstacle_Area,Covered_Area);
                alpop_Cost(2) = Energy_consumption(alpop,rc);
                alpop_Cost(3) = Late_Func(alpop,rc);
                
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
        pop(i).Velocity(j,1:2) = pop(i).Velocity(j,1:2) * wdamp;
    end

    disp(['Iteration ' num2str(it) ':' num2str(size(non_dom_pop,1)) ' Non-dominated solutions']);
    clear h j k phi alpop alpop2 alpop_Cost E E_ave;   
end
%% plot
% %alpop=pop(4).Position;
% alpop=BestSol.Position;
% Cost_Function = Cost_Func(alpop,weight,stat,Obstacle_Area,Covered_Area);
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
% %clear i x_fill y_fill theta obs_col obs_row;
% 
% axis equal;
% xlim([0 size(Obstacle_Area,1)]);
% ylim([0 size(Obstacle_Area,2)]);
% title(['Weighted Coverage ratio:' num2str(BestSol.Cost*100) '% ' ]);
% grid on;
% drawnow;
%{
if save_fig==1
    saveas(gcf, ['Iteration' num2str(it) '.pdf' ]); % Lưu dưới dạng PDF
end
%}
%%
%save(name)
%end
figure;
hold on;
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
