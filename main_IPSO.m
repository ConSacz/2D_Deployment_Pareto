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
MaxIt = 200;              % Maximum Number of Iterations
N = 25;
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
BestSol.Position=[];	% Best Solution
BestSol.Cost=1;	% Best solution function value
pop=repmat(empty_individual,nPop,1); % pop includes nPop solutions
BestCost=zeros(MaxIt,1); % Best solution every iteration

clear empty_individual;

for k = 1:nPop
    alpop = unifrnd(sink(1)-rc,sink(2)+rc,[N 2]);
    % gen sink node
    alpop (1,:)= sink;
    pop(k).Position=alpop;
    pop(k).Cost=Cost_Func(pop(k).Position,weight,stat,Obstacle_Area,Covered_Area);
    pop(k).Velocity=zeros(N,3);
    % check for best Solution
    if pop(k).Cost < BestSol.Cost
	    BestSol.Cost=pop(k).Cost;
	    BestSol.Position=pop(k).Position;
    end 
end
clear k alpop sink;
%%     GA Main Loop
for it = 1:MaxIt
    for i=1:nPop
        for j=1:N
            alpop=pop(i).Position;
            k=randi([1,N]);
            alpop_v=w*pop(i).Velocity(j,1:2) + c1*rand()*(pop(i).Position(k,1:2) - pop(i).Position(j,1:2)) + c2*rand()*(BestSol.Position(j,1:2) - pop(i).Position(j,1:2));
            alpop(j,1:2) = pop(i).Position(j,1:2) + alpop_v;
            alpop(:,1:2) = min(max(alpop(:,1:2), 1),size(Obstacle_Area,1)-1);
            if Connectivity_graph(Graph(alpop(:,1:2),rc),[]) == 1
                %[pop_cov,~] = Cov_Func(pop(i).Position,rs,theta0,Obstacle_Area,Covered_Area);
                alpop_cov  = Cost_Func(alpop,weight,stat,Obstacle_Area,Covered_Area);
                if alpop_cov <= pop(i).Cost
                    pop(i).Position = alpop;
                    pop(i).Cost = alpop_cov;
                    pop(i).Velocity(j,1:2) = alpop_v;
                end
            end
        end
        pop(i).Velocity(j,1:2) = pop(i).Velocity(j,1:2) * wdamp;
    end


    clear alpop_cov alpop alpop_v alpop2;
    for i=1:nPop
        if pop(i).Cost < BestSol.Cost
            BestSol.Cost=pop(i).Cost;
            BestSol.Position=pop(i).Position;
        end 
    end
    BestCost(it)=BestSol.Cost;
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
    clear h j k phi alpop alpop2 alpop_cov E E_ave;   
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

