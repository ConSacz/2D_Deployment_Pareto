%% Directional Sensor
%% DEPLOYMENT

%%
clc;
clear;

%% 
close all;
%% Network parameter
% for N=30:10:50
% for trial=1:20

% Monitor area
Covered_Area = zeros(100,100);
%Obstacle_Area = gen_random_distribution_area();
Obstacle_Area = ones(100,100);
%Obstacle_Area=genarea();
%Obstacle_Area = gen_target_area(1000);

% nodes info
MaxIt = 1000;              % Maximum Number of Iterations
a = 1;                    % Acceleration Coefficient Upper Bound
N = 60;
nPop=25;
rc = 10; 
rs=ones(1,N)*10;
stat(1,:)=rs;
stat(2,1)=rc;
weight=[0.5 0.5 0.5];
%name=['./data/hetero_target/',num2str(N),' nodes/case4/hetero_',num2str(trial),'.mat'];

sink=[50 50];
Scout_bee=nPop;
Onlooker_bee=nPop;

% bee parameter
a=1;


%% Init first pop
empty_individual.Position=[]; % Empty solution
empty_individual.Cost=[];     % Empty cost function of that solution
BestSol.Position=[];	% Best Solution
BestSol.Cost=1;	% Best solution function value
pop=repmat(empty_individual,nPop,1); % pop includes nPop solutions
BestCost=zeros(MaxIt,1); % Best solution every iteration
L=zeros(1,nPop);

clear empty_individual;

for k = 1:nPop
    alpop = unifrnd(sink(1)-rc/2,sink(2)+rc/2,[N 2]);
    % gen sink node
    alpop (1,:)= sink;
    pop(k).Position=alpop;
    pop(k).Cost=Cost_Func(pop(k).Position,weight,stat,Obstacle_Area,Covered_Area);

    % check for best Solution
    if pop(k).Cost<BestSol.Cost
	    BestSol.Cost=pop(k).Cost;
	    BestSol.Position=pop(k).Position;
    end 
end
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
            %[pop_cov,~] = Cov_Func(pop(i).Position,rs,theta0,Obstacle_Area,Covered_Area);
            alpop_cov  = Cost_Func(alpop,weight,stat,Obstacle_Area,Covered_Area);
            if alpop_cov <= pop(i).Cost
                pop(i).Position=alpop;
                pop(i).Cost=alpop_cov;
            else
                L(i)=L(i)+1;
            end
            break;
        end
    end

    clear alpop_cov alpop alpop2 k phi;
    %% ranking pop
    E=zeros(1,nPop);
    for i=1:nPop
        E(i)=pop(i).Cost;
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
        % change the position and angle of each node
        for k=1:N
            alpop=pop(i).Position;
            h=randi([1,N]);
            phi=a*unifrnd(-1, +1, [1 2])*(1-L(i)/MaxIt)^2;
            alpop(k,1:2)  = pop(i).Position(k,1:2) + phi.*(pop(i).Position(k,1:2)-pop(i).Position(h,1:2));
            alpop(:,1:2) = min(max(alpop(:,1:2), 1),size(Obstacle_Area,1)-1);
            if Connectivity_graph(Graph(alpop(:,1:2),rc),[]) == 1
                %[pop_cov,~] = Cov_Func(pop(i).Position,rs,theta0,Obstacle_Area,Covered_Area);
                alpop_cov  = Cost_Func(alpop,weight,stat,Obstacle_Area,Covered_Area);
                if alpop_cov <= pop(i).Cost
                    pop(i).Position=alpop;
                    pop(i).Cost=alpop_cov;
%                 else
%                     L(i)=L(i)+1;
                end
            end
        end
    end
    clear alpop_cov alpop alpop2 k h phi;
    for i=1:nPop
        if pop(i).Cost<BestSol.Cost
            BestSol.Cost=pop(i).Cost;
            BestSol.Position=pop(i).Position;
        end 
    end
    BestCost(it)=BestSol.Cost;
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
    clear h j k phi alpop alpop2 alpop_cov E E_ave;
    if BestCost(it)==1
        break;
    end    
end
%% plot
G=Graph(BestSol.Position,rc);
clf();
hold on;
for i = 1:1:numel(G.Edges.EndNodes)/2
    plot([BestSol.Position(G.Edges.EndNodes(i,1),1),BestSol.Position(G.Edges.EndNodes(i,2),1)],[BestSol.Position(G.Edges.EndNodes(i,1),2),BestSol.Position(G.Edges.EndNodes(i,2),2)],'Color','blue','linewidth',1);
end
for i = 1:N
    plot (BestSol.Position(i,1) , BestSol.Position(i,2),'ro');
    %viscircles ([BestSol.Position(i,1) BestSol.Position(i,2)],rs(i),'Color', 'k');
    %text (BestSol.Position(i,1) , BestSol.Position(i,2), num2str(i),'FontSize',15,'Color','red');
end

%plot(obs_row, obs_col,'.', 'MarkerSize', 20, 'Color', 'red');
xlim([0 100]);
ylim([0 100]);
title(['Coverage Ratio: ', num2str(100-BestSol.Cost*100),'%']);
grid on;


%%
% save(name)
% end
% end
