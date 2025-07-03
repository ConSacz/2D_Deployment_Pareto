function [Cost] = Sphere_MO(pop, Rs, Area)
a=0.5;
b=0.5;
Cov=Sphere(pop,Rs,Area);
Life=Life_Time(Graph(pop,Rs));
%{
if(Cov>0.25)
    Cost= a*Cov + b*Life;
else
    Cost= a*Cov;
end
%}
Cost= a*Cov + b*Life;