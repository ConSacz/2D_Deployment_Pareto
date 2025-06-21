function [Cost] = CostFunction(pop,stat,Obstacle_Area,Covered_Area)

rs=stat(1,:);
rc=stat(2,1);

G = Graph(pop,rc);

Coverage = Cov_Func_v2(pop,rs,Obstacle_Area,Covered_Area);

Energy = Energy_consumption(G);

LifeTime = Life_Time(G);


Cost= [Coverage  Energy  LifeTime]';