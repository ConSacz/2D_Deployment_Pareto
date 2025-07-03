clc;
clear;
Rs=10;
Rc=10;
rawData.Cost=[];
minIt=100;
maxIt=5300;
nPop=200;
allpop=(maxIt/100-minIt/100+1)*nPop;
Data=repmat(rawData,allpop,1);
clear rawData;
%% load data
%load ('ran_pop.mat')
for j=minIt:100:maxIt
    name =[ num2str(j) '.mat'];
    load (name)
    for i =1:nPop
        G=Graph(pop(i).Position,Rc);
        %Data((j-minIt)/100*nPop+i).Cost(1)=Sphere(pop(i).Position,Rs,[100 100]) ;
        %Data((j-minIt)/100*nPop+i).Cost(2)=Energy_consumption(G);
        Data((j-minIt)/100*nPop+i).Cost(3)=Life_Time(G);
    end
    disp(['Iteration ' num2str(j) ' done ']); 
end
%% pareto front rank 1 la duong chua cac loi giai dominate
a1=[1];                 %ma tran chua cac vi tri cua loi giai dominate
for i=2:allpop          %chay giua cac diem trong Data
    for j=1:numel(a1)
        % if a1 dominates pop(i)
        if checkDomination(Data(a1(j)).Cost,Data(i).Cost) == 1
            count=0;
            break;
        % if pop(i) dominated a1
        elseif checkDomination(Data(a1(j)).Cost,Data(i).Cost) == -1
            count=0;
            a1(j)=i;
        % nondominated
        else
            count=1;
        end
    end
    if count ==1
        a1=[a1 i];
    end
    a1=unique(a1);
end

%% Decision maker weighted metrics
% decision=[];
% a1N=numel(a1);
% a=3;
% w1=[0.5 0.5];
% zs=[0 0];
% wmetric=zeros(1,a1N);
% for i=1:a1N
%     wmetric(i)=( w1(1)*(Data(a1(i)).Cost1-zs(1))^a + w1(2)*(Data(a1(i)).Cost2-zs(2))^a )^(1/a);
% end
% [val1, minidx]=min(wmetric);
% [val2, maxidx]=max(wmetric);
% decision=[decision a1(minidx)];
% decision=[decision a1(maxidx)];

%% Decision maker Tchebychef
% Tchebychef=zeros(2,a1N);
% w2=[0.5 0.5];
% for i=1:a1N
%     Tchebychef(1,i)=w2(1)*(Data(a1(i)).Cost1-zs(1));
%     Tchebychef(2,i)=w2(2)*(Data(a1(i)).Cost2-zs(2));
% end
% [val1, minidx2]=min(Tchebychef(1,:));
% [val2, minidx3]=min(Tchebychef(2,:));
% minidx2=a1(minidx2);
% minidx3=a1(minidx3);
%decision=[decision a1(minidx)];
%decision=[decision a1(maxidx)];
%%
%{
PRT1=zeros(allpop,1);
PRT2=zeros(allpop,1);

for i =1:allpop
    for j =1:allpop
         if ((Data(i).Cost1==Data(j).Cost1)&&(Data(i).Cost2>Data(j).Cost2)) || ((Data(i).Cost1>Data(j).Cost1)&&(Data(i).Cost2==Data(j).Cost2)) || ((Data(i).Cost1>Data(j).Cost1)&&(Data(i).Cost2>Data(j).Cost2))
         %if (Data(i).Cost1<=Data(j).Cost1)&&(Data(i).Cost2<=Data(j).Cost2&&i~=j)    
             PRT1(i)=PRT1(i)+1;
             PRT2(j)=PRT2(j)+1;
         end
    end
end
%}
%% plot PF
figure;
% xlim([0 1])
% ylim([0 1])
%hold on;
for i=1:numel(a1)
     %if a1(i)== 4967
     %   plot (Data(a1(i)).Cost1 , (Data(a1(i)).Cost2),'ro','Color','r');
     %elseif a1(i)== decision(2)
     %    plot (Data(a1(i)).Cost1 , (Data(a1(i)).Cost2),'ro','Color','b');
     %else
        plot3 (Data(a1(i)).Cost(1),Data(a1(i)).Cost(2),Data(a1(i)).Cost(3), 'ro','Color','g');
     %end
     
     %text (Data(a1(i)).Cost1 , Data(a1(i)).Cost2, num2str(a1(i)),"FontSize",6)
     hold on;
end
xlabel('Non-coverage cost function') 
ylabel('Energy cost function')
zlabel('Life time cost function')
axis equal;
xlim([0 1]);
ylim([0 1]);
xticks(0:0.1:1);
yticks(0:0.1:1);

%% plot all
f = figure;
%f.Position = [100 100 1000 600];
%axis([0 1 0 3])
for i=1:allpop
    data(i,1)=Data(i).Cost1;
    data(i,2)=Data(i).Cost2;
    %ran_data(i,1)=ran_Data(i).Cost1;
    %ran_data(i,2)=ran_Data(i).Cost2;
end

p1=plot (data(:,1) , data(:,2),'o','Color','g');
hold on;
%p2=plot (ran_data(:,1) , ran_data(:,2),'o','Color','c');
for i=1:numel(a1)
    r1_data(i,1)=Data(a1(i)).Cost1;
    r1_data(i,2)=Data(a1(i)).Cost2;
end
plot (r1_data(:,1), r1_data(:,2),'o','Color','r');
%plot (ran_Data(a2(i)).Cost1 , ran_Data(a2(i)).Cost2,'o','Color','b');

legend('pareto set','nondominant solutions')
xlabel('Non-coverage cost function') 
ylabel('Energy cost function') 
ax = gca;
exportgraphics(ax,"rank1.pdf","Resolution",300)
%% plot decision
%get the pop of decision
for j=minIt:100:maxIt
    load (num2str(j)+".mat")
    for i =1:nPop
        if ((j-minIt)/100*nPop+i)==3839
            pop1=pop(i).Position;
        end
    end
    %disp(['Iteration ' num2str(j) ' done ']); 
end

% Graph
adj_pop=zeros(numel(pop1)/2);
for i=1:(numel(pop1)/2)
    for j=1:(numel(pop1)/2)
        dist = sqrt((pop1(i*2)-pop1(j*2))^2+(pop1(i*2-1)-pop1(j*2-1))^2);
        if dist <= Rs && dist ~=0
            adj_pop(i,j)=dist;
        end
    end
end
G=graph(adj_pop);
clear i j dist adj_pop ;

% plot nodes and sensing range
BestSol.Position=pop1;
figure;
hold on;

for i = 1:numel(BestSol.Position)/2   
    %hold on;
    %text (BestSol.Position(1,i) , BestSol.Position(1,i+1), num2str(i/2+0.5),'FontSize',25);
    %hold on;
    %viscircles ([BestSol.Position(1,i*2) BestSol.Position(1,i*2-1)],Rs,'Color', 'k','linewidth',0.1);
    p = nsidedpoly(1000, 'Center', [BestSol.Position(1,i*2) BestSol.Position(1,i*2-1)], 'Radius', Rs);
    plot(p, 'FaceColor', 'c')
    
end

for i = 1:numel(BestSol.Position)/2  
    plot (BestSol.Position(1,i*2) , BestSol.Position(1,i*2-1),'o','Color','r');
end
%text (BestSol.Position(1,2) , BestSol.Position(1,3), num2str(1/2+0.5),'FontSize',25);
grid on;



title('Coverage Ratio = 92.21%')
%plot paths
%{
for i=2:(numel(pop1)/2)
    path = shortestpath(G,i,1);
    for member = 1:(numel(path)-1)
        p1 = [pop1(1,path(member)*2)   pop1(1,path(member)*2-1)  ];                         % First Point
        p2 = [pop1(1,path(member+1)*2) pop1(1,path(member+1)*2-1)];                         % Second Point
        dp = p2-p1;                                                                         % Difference
        quiver(p1(1),p1(2),dp(1),dp(2),0,'Color','b','LineWidth',0.7,'MaxHeadSize',1.2)
    end
end
%}
axis equal;
xlim([0 100]);
ylim([0 100]);
xticks(0:10:100);
yticks(0:10:100);
%drawnow;
ax = gca;
exportgraphics(ax,"plot_decision_cr2.pdf","Resolution",300)

%%
figure
axes('Box','on');
%tpos = tightPosition(axes);
%annotation("rectangle",tpos,LineWidth=0.3)
for i=1:allpop
    data(i,1)=Data(i).Cost1;
    data(i,2)=Data(i).Cost2;
    ran_data(i,1)=ran_Data(i).Cost1;
    ran_data(i,2)=ran_Data(i).Cost2;
end
hold on;

p=plot (data(:,1) , data(:,2),'o','Color','g');
%p2=plot (data(1001:2000,1) , data(1001:2000,2),'o','Color','r');
%p3=plot (data(2001:3000,1) , data(2001:3000,2),'o','Color','b');
%p4=plot (data(3001:4000,1) , data(3001:4000,2),'o','Color','k');
p5=plot (data(4001:5000,1) , data(4001:5000,2),'o','Color','b');
%ran_p=plot (ran_data(:,1) , ran_data(:,2),'o','Color','c');
xlim([0 1])
ylim([0 1.8])
legend('weighted-sum based ABC pareto set','traditional ABC pareto set')
%legend('random pareto set')
xlabel('Non-coverage cost function') 
ylabel('Energy cost function')
ax = gca;
exportgraphics(ax,"pareto sets.pdf","Resolution",300)