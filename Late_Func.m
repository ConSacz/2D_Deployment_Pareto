function Latency = Late_Func(pop,rc)
% Latency is a COST FUNCTION
% this function calculate the max value of latency 

%% Network parameter
N = size(pop,1);
c = 0.001;

%% Generate network Graph
adj_pop=zeros(N);
for i=1:(N)
    for j=1:(N)
        dist = sqrt((alpop(i,1)-alpop(j,1))^2+(alpop(i,2)-alpop(j,2))^2);
        if dist <= rc && dist ~=0
            adj_pop(i,j)=dist;
        end
    end
end
G=graph(adj_pop);

%% Calculate each node latency when transfering data to sink
Late=zeros(N,1);
for j=2:N
    path = shortestpath(G,1,j);
    
    Late(j) = c * numel(path);
end

%% Get the maximum latency
Latency = max(Late);