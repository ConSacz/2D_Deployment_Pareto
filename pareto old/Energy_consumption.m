function E_normalized = Energy_consumption(G)
% Energy Consumption is a COST FUNCTION
% This function calculate total energy consumption of a network in a
% operation round
% min = 0.28, max = 5.49 in 60-node 10-10 case
%% parameter
N=numnodes(G);
b=0.1;  % nJ/bit/m^a environment constant
a=2;    % exponent path loss
EM=13;  % nJ/bit    maintain/process energy
ET=20;  % nJ/bit    transmit energy
ER=2;   % nJ/bit    receive energy
maxBat=10000;

% normalize
minval=0.2;
maxval=6;

%% Calculate each node energy consumption using route as shortest path
Bat=zeros(N,1);
for j=2:N
    path = shortestpath(G,1,j);
    for i=1:size(path,2)
        if i==1
           %Bat(path(i)) = Bat(path(i))-(EM+ET) ;
           %Bat(path(i))=0;
           continue;
        elseif i==size(path,2)
           dt = G.Edges.Weight(findedge(G,path(i),path(i-1)));
           Bat(path(i)) = Bat(path(i))+(2*EM+ET+b*dt^a) ; 
        else
           dt = G.Edges.Weight(findedge(G,path(i),path(i-1)));
           dr = G.Edges.Weight(findedge(G,path(i),path(i+1)));
           Bat(path(i)) = Bat(path(i))+(EM+ER+ET+b*dt^a+b*dr^a) ; 
        end
    end
end

%% Sum up 
E_total=sum(Bat)/maxBat;
% normalize 0-inf to 0-1
E_normalized = round((E_total-minval)/(maxval-minval),2);
%E_normalized = round(E_total/(1+E_total),3);