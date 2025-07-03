function coverage = Sphere(pop,Rs,Area)
%x is a 2dim matrix holding position of nodes
%Radius is the sensing rad of nodes
%area is a 2dim matrix, stands for the measuring area

count=0;
pts=[100,100];     %distribute pts points throughout the map
%----------------- create the point map
%pointspos=[0:(area(1)/(pts/2-1)):area(1);0:(area(2)/(pts/2-1)):area(2)];
%-----------------
%coveragarea= zeros(pts(1),pts(2));
for i=0:(pts(1))
    for k=0:(pts(2))
        for j=1:2:numel(pop)
            dist = sqrt((i*Area(1)/pts(1)-pop(1,j))^2+(k*Area(2)/pts(2)-pop(1,j+1))^2);
            if dist< Rs || dist== Rs
                %coveragarea(i,k)=1;
                count=count+1;
                break   
            end
        end
    end    
end
coverage=1-((count)/((pts(1)+1)*(pts(2)+1)));


