function [Cost] = Cost_Func(pop,weight,stat,Obstacle_Area,Covered_Area)

rs=stat(1,:);
rc=stat(2,1);
a=weight(1);
b=weight(2);
c=weight(3);

G = Graph(pop,rc);

if a ~= 0
    Coverage = a * Cov_Func_v2(pop,rs,Obstacle_Area,Covered_Area);
else 
    Coverage = 0; 
end
if b ~= 0
    Energy = b * Energy_consumption(G);
else 
    Energy = 0; 
end
if c ~= 0
    LifeTime = c * Life_Time(G);
else 
    LifeTime = 0; 
end

Cost= Coverage + Energy + LifeTime;