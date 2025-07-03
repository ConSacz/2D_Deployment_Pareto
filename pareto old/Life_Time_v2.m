function Life_Ratio = Life_Time_v2(pop,Rs)
N=size(pop,2)/2;
Life_Ratio=0;
adj_pop=zeros(size(pop,2)/2);
for i=1:(size(pop,2)/2)
    for j=1:(size(pop,2)/2)
        dist = sqrt((pop(i*2)-pop(j*2))^2+(pop(i*2-1)-pop(j*2-1))^2);
        if dist <= Rs && dist ~=0
            adj_pop(i,j)=1;
        end
    end
end
G=graph(adj_pop);
clear i j dist adj_pop ;
%%

nb=neighbors(G,1);
count=zeros(1,numel(nb));
for i=1:N
    path=shortestpath(G,i,1);
    for j=1:numel(nb)
        if ismember(nb(j),path)
            count(j)=count(j)+1;
        end
    end
end
clear i j;
max_node=max(count);
for i=1:numel(nb)
    Life_Ratio = Life_Ratio+count(i)/max_node;
end
Life_Ratio=1/round(Life_Ratio,6);
%%