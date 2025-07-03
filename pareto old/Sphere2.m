function coverage = Sphere2(pop,Rs,Area,Area1)
pop = reshape(pop,[numel(pop)/2,2]);
count1=0;
count2=0;
pts=[100,100];     %distribute pts points throughout the map
%----------------- create the point map
%pointspos=[0:(area(1)/(pts/2-1)):area(1);0:(area(2)/(pts/2-1)):area(2)];
%-----------------
%coverarea= zeros(pts(1),pts(2));
Area=[100,100];
coverarea=Area1;
for i=1:(pts(1))
    for k=1:(pts(2))
        for j=1:size(pop,1)
            dist = sqrt((i*Area(1)/pts(1)-pop(j,1))^2+(k*Area(2)/pts(2)-pop(j,2))^2);
            if (dist< Rs || dist== Rs) && coverarea(i,k)== 0
                coverarea(i,k)=1;
                count1=count1+1;
                break   
            end
            if (dist< Rs || dist== Rs) && coverarea(i,k)== -1
                coverarea(i,k)=-2;
                count2=count2+1;
                break   
            end
        end
    end    
end
coverage=1/((count1-count2)/(pts(1)*pts(2)));


