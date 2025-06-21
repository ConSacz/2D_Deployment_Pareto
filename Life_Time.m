function lifetime_normalized = Life_Time(G)
% Life time is a FITNESS FUNCTION
% This function calculate life time of a network
% min = 0.05, max = 2.9 in 60-node 10-10 case
%% parameter
N=numnodes(G);
b=0.1;  % nJ/bit/m^a environment constant
a=2;    % exponent path loss
EM=13;  % nJ/bit    maintain/process energy
ET=20;  % nJ/bit    transmit energy
ER=2;   % nJ/bit    receive energy
maxBat=10000;

% normalize
minval = 0;
maxval = 4;

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
lifetime= maxBat/max(Bat);
% normalized 0-inf to 1-0
%lifetime_normalized = round((10/lifetime-minval)/(maxval-minval),2);
lifetime_normalized = round(1/lifetime,3);
%lifetime_normalized = round(1/(1+lifetime),2);