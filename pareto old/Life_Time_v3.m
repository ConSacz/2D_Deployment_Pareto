function Life_Ratio = Life_Time_v3(pop,Rs)
%Life time COST FUNCTION need to min
N=size(pop,2)/2;
EM=13;
ET=20;
ER=2;
maxBat=1000;

adj_pop=zeros(size(pop,2)/2);
for i=1:(size(pop,2)/2)
    for j=1:(size(pop,2)/2)
        dist = sqrt((pop(i*2)-pop(j*2))^2+(pop(i*2-1)-pop(j*2-1))^2);
        if dist <= Rs && dist ~=0
            adj_pop(i,j)=dist;
        end
    end
end
G=graph(adj_pop);
clear i j dist adj_pop ;
%%
Bat=zeros(N,1);
for j=2:N
    path = shortestpath(G,1,j);
    for i=1:size(path,2)
        if i==1
           %Bat(path(i)) = Bat(path(i))-(EM+ET) ;
           Bat(path(i))=0;
        elseif i==size(path,2)
           Bat(path(i)) = Bat(path(i))+(EM+ER) ; 
        else
           Bat(path(i)) = Bat(path(i))+(EM+ER+ET) ; 
        end
    end
end
%%
Life_Ratio=max(Bat)/maxBat;