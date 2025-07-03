clc;
clear;

N=60;
nPop=10000;

VarMin = 0;         % Decision Variables Lower Bound

VarMaxx = 100;         % Decision Variables Upper Bound
VarMaxy = 100;
Area = [VarMaxx VarMaxy]; %Area [x1 x2 .. xn y1 .. yn]

nVar = N * 2;             % Variables' dimensions

VarSize = [1 nVar];   % Decision Variables Matrix Size

Rs = 10; %Transmission range = 10 meter
Rc = 10; %Communication range = 10 meter

rawData.Cost1=[];
rawData.Cost2=[];
Data=repmat(rawData,nPop,1);
clear rawData
%% Random Initialization and Waggle dance
% pop gen
pop=zeros(nPop,nVar);
for k = 1:nPop
    pos = zeros(1,2*N);
    pos (1)= VarMaxx/2;
    pos (2)= VarMaxy/2;
    for i=2:N
        Rcom= Rc;
        j=randi(i-1,1);
        pos (i*2-1)= pos(j*2-1)+2*Rcom*rand-Rcom;
        if (rand>0.5)
            pos (i*2)= sqrt(Rcom^2-(pos(i*2-1)-pos(j*2-1))^2)+pos(j*2);
        else
            pos (i*2)= -sqrt(Rcom^2-(pos(i*2-1)-pos(j*2-1))^2)+pos(j*2);
        end
    end
    pop(k,:) =  pos; %unifrnd(VarMin, VarMax, VarSize);
    % update cost
    Data(k).Cost1=Sphere(pop(k,:),Rs,Area) ;
    Data(k).Cost2=Life_Time_v3(pop(k,:),Rc);
end 
%% sort rank 1
a1=[1];                 %ma tran chua cac vi tri cua loi giai dominate
for i=2:nPop          %chay giua cac diem trong Data
    for j=1:numel(a1)
        %neu ca 2 CF deu lom hon hoac bang
        if (Data(i).Cost1>=Data(a1(j)).Cost1)&&(Data(i).Cost2>=Data(a1(j)).Cost2)
            count=0;
            break;
        %neu 1 cai xin hon han, 1 cai xin hon hoac bang thi thay the cai do 
        elseif ((Data(i).Cost1<=Data(a1(j)).Cost1)&&(Data(i).Cost2<Data(a1(j)).Cost2)) || ((Data(i).Cost1<Data(a1(j)).Cost1)&&(Data(i).Cost2<=Data(a1(j)).Cost2))
            count=0;
            a1(j)=i;
        %ko xin hon cung ko lom hon
        else
            count=1;
        end
    end
    if count ==1
        a1=[a1 i];
    end
    a1=unique(a1);
end
%% plot rank 1
%{
figure;
for i=1:numel(a1)
     plot (Data(a1(i)).Cost1 , (Data(a1(i)).Cost2),'ro','Color','r');
     %text (Data(a1(i)).Cost1 , 1/Data(a1(i)).Cost2, num2str(a1(i)),"FontSize",6)
     text (Data(a1(i)).Cost1 , Data(a1(i)).Cost2, num2str(a1(i)),"FontSize",6)
     hold on;
end
%}
%% plot all
figure;
axis([0 1 0 3])
for i=1:nPop
    if ismember(i,a1)
        plot (Data(i).Cost1 , Data(i).Cost2,'ro','Color','r');
    else
        plot (Data(i).Cost1 , Data(i).Cost2,'ro','Color','b');
    end
    hold on;
end
%%
a2=a1;
ran_pop=pop;
ran_Data=Data;
save("ran_pop.mat","a2","ran_pop","ran_Data")