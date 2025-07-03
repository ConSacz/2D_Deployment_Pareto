function Result = Connectivity_v2(pop,Rc)
% result=0 so connected
numberNodes = numel(pop)/2;
count=zeros(numberNodes,1);
sec=0;
Result=1;
for i= 1:numberNodes
    for j= 1:numberNodes
        if (((pop(2*i-1)-pop(2*j-1))^2+(pop(2*i)-pop(2*j))^2)<=(Rc)^2)
            if count(i)==0 && count(j)==0
                sec=sec+1;
                count(j)=sec;
                count(i)=count(j);
            elseif count (i) == 0 && count(j)~=0
                count(i)=count(j);   
            elseif count (i) ~= 0 && count(j)==0
                count(j)=count(i);    
            elseif count (i) ~= 0 && count(j)~=0
                count(i)= min(count(i),count(j));
                count(j)=count(i);
            end
        end
    end
end
if count == ones(size(pop,1),1)
    Result=0;
end