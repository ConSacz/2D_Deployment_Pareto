function G=Graph(pop,rc)
N=numel(pop)/2;
adj_pop=zeros(size(pop,2)/2);
for i=1:(size(pop,2)/2)
    for j=1:(size(pop,2)/2)
        dist = sqrt((pop(i*2)-pop(j*2))^2+(pop(i*2-1)-pop(j*2-1))^2);
        if dist <= rc && dist ~=0
            adj_pop(i,j)=1;
        end
    end
end
G=graph(adj_pop);
