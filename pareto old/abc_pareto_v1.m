clc;
clear;
N = 60 ; %Number of sensor nodes

VarMin = 0;         % Decision Variables Lower Bound

VarMaxx = 100;         % Decision Variables Upper Bound
VarMaxy = 100;
Area = [VarMaxx VarMaxy]; %Area [x1 x2 .. xn y1 .. yn]

nVar = N * 2;             % Variables' dimensions

VarSize = [1 nVar];   % Decision Variables Matrix Size

Rs = 10; %Transmission range = 10 meter

%% ABC Settings

MaxIt = 2000;              % Maximum Number of Iterations

nPop = 200;               % Population Size (Colony Size)

nOnlooker = nPop;         % Number of Onlooker Bees

Limit = nPop * nVar; % Abandonment Limit Parameter (Trial Limit)

a = 1;                    % Acceleration Coefficient Upper Bound

Rc =  Rs;
%% Random Initialization and Waggle dance
% Empty Bee Structure
empty_bee.Position = [];
empty_bee.Cost = [];

% Initialize Population Array
pop = repmat(empty_bee, nPop, 1);

% Initialize Best Solution Ever Found
BestSol.Cost = Inf;
BestSol.Position=zeros(1,2*N);
% Create Initial Population
for k = 1:nPop
    pos = zeros(1,2*N);
    pos (1)= VarMaxx/2;
    pos (2)= VarMaxy/2;
    for i=2:N
        Rcom= Rc*rand;
        j=randi(i-1,1);
        pos (i*2-1)= pos(j*2-1)+2*Rcom*rand-Rcom;
        if (rand>0.5)
            pos (i*2)= sqrt(Rcom^2-(pos(i*2-1)-pos(j*2-1))^2)+pos(j*2);
        else
            pos (i*2)= -sqrt(Rcom^2-(pos(i*2-1)-pos(j*2-1))^2)+pos(j*2);
        end
    end
    pop(k).Position =  pos; %unifrnd(VarMin, VarMax, VarSize);
     CostFunction1 = Sphere_MO(pop(k).Position, Rs, Area);
     pop(k).Cost = CostFunction1;
     if pop(k).Cost <= BestSol.Cost
         BestSol = pop(k);
     end
end 
% Abandonment Counter
 C = zeros(nPop, 1);
 
 % Array to Hold Best Cost Values
 BestCost = zeros(MaxIt, 1);
 BestFit = zeros(MaxIt,1);
 %% Local search
for it = 1:MaxIt
    for i = 1:nPop
        while (1)
            % Choose k randomly, not equal to i
            K = [1:i-1 i+1:nPop];
            k = K(randi([1 numel(K)])); %Randomly neibourhood bee
        
            % Define Acceleration Coeff.
            phi = a*unifrnd(-1, +1);
            w1 = randi([1 N]);
            w1=w1*2-1;
            % New Bee Position
            w1_new=w1+1;
            newbee.Position = pop(i).Position;
            newbee.Position(w1) = pop(i).Position(w1)+phi.*(pop(i).Position(w1)-pop(k).Position(w1));
            newbee.Position(w1_new) = pop(i).Position(w1_new)+phi.*(pop(i).Position(w1_new)-pop(k).Position(w1_new));
            % Apply Bounds
            newbee.Position = max(newbee.Position, VarMin);
            newbee.Position(1,1:2:2*N-1) = min(newbee.Position(1,1:2:2*N-1), VarMaxx);
            newbee.Position(1,2:2:2*N) = min(newbee.Position(1,2:2:2*N), VarMaxy);
            check1 = Connectivity(newbee.Position, Rc);
            if check1 == 0
                newbee.Position = pop(i).Position;
            else
                break;
            end
        end
    end
    
    % Calculate Fitness Values and Selection Probabilities
    F = zeros(nPop, 1);
    MeanCost = mean([pop.Cost]);
    for i = 1:nPop
        %-------------------------------------------------------------------
        F(i) = exp(-pop(i).Cost/MeanCost); % Convert Cost to Fitness
        %-------------------------------------------------------------------
    end
    P = F/sum(F);
    
    % Onlooker Bees
    for m = 1:nOnlooker % run m times but no use m as a matrix counter
        
        % Select Source Site
        %-------------------------------------------------------------------
        i = RouletteWheelSelection(P); %not run i from 1 to pop but randomly select i using roullette/greedy
        %-------------------------------------------------------------------
        % Choose k randomly, not equal to i
        while (1)
            K = [1:i-1 i+1:nPop];
            k = K(randi([1 numel(K)]));
        
        % Define Acceleration Coeff.
            phi = a*unifrnd(-1, +1);
            w2 = randi([1 N]);
            w2 = w2*2-1;
            w2_new = w2+1;
        % New Bee Position
            newbee.Position = pop(i).Position; 
            newbee.Position(w2) = pop(i).Position(w2)+phi.*(pop(i).Position(w2)-pop(k).Position(w2));
            newbee.Position(w2_new) = pop(i).Position(w2_new)+phi.*(pop(i).Position(w2_new)-pop(k).Position(w2_new));
        % Apply Bounds
            newbee.Position = max(newbee.Position, VarMin);
            newbee.Position(1,1:2:2*N-1) = min(newbee.Position(1,1:2:2*N-1), VarMaxx);
            newbee.Position(1,2:2:2*N) = min(newbee.Position(1,2:2:2*N), VarMaxy);
            check1 = Connectivity(newbee.Position, Rc);
            if check1 == 0
                newbee.Position = pop(i).Position;
            else
                break;
            end
        end

        % Evaluation
        CostFunction2 = Sphere_MO(newbee.Position, Rs, Area);
        newbee.Cost = CostFunction2; 
        
        % Comparision
        if newbee.Cost <= pop(i).Cost
            pop(i) = newbee;
        else
            C(i) = C(i) + 1;
        end
        
    end
    %-----------------------------------------------------------------------
    %% Global Search
    % Scout Bees
    for k = 1:nPop
        if C(k) >= Limit
                pos = zeros(1,2*N);
                pos (1)= VarMaxx/2;
                pos (2)= VarMaxy/2;
                for i=2:N
                    Rcom= Rc*rand;
                    j=randi(i-1,1);
                    pos (i*2-1)= pos(j*2-1)+2*Rcom*rand-Rcom;
                    if (rand>0.5)
                        pos (i*2)= sqrt(Rcom^2-(pos(i*2-1)-pos(j*2-1))^2)+pos(j*2);
                    else
                        pos (i*2)= -sqrt(Rcom^2-(pos(i*2-1)-pos(j*2-1))^2)+pos(j*2);
                    end
                end
                
                pop(k).Position =  pos; %unifrnd(VarMin, VarMax, VarSize);
                CostFunction3 = Sphere_MO(pop(k).Position, Rs, Area);
                pop(k).Cost = CostFunction3;
                C(k) = 0;
        end
    end
    %-----------------------------------------------------------------------delete
    %exhausted source out of loop
    
    % Update Best Solution Ever Found
    for i = 1:nPop
        if pop(i).Cost <= BestSol.Cost
            BestSol = pop(i);
        end
    end 

    if mod(it,2000)==0
        %pop1=pop(1:200);
        %pop2=pop(201:400);
        %save(num2str(it*100+500)+".mat",'pop1');
        %save(num2str(it*100+800)+".mat",'pop2');
        %save(num2str(it/2+8000)+".mat",'pop');
        save("5200.mat",'pop');
        %clear pop1 pop2;
    end
    % Store Best Cost Ever Found
    BestCost(it) = BestSol.Cost;
    BestPoisition = BestSol.Position;
    BestFit(it) = exp(-BestSol.Cost/mean([pop.Cost]));
    % Display Iteration Information
    disp(['Iteration ' num2str(it) ': Best ratio coverage = ' num2str(BestCost(it))]);   
end
%for i = 1:nPop
    %disp([ 'x' : num2str(pop(i).Position)]);
%end
%% Results
%{
figure;
semilogy(BestCost,'LineWidth',2);
% plot(BestCost,'LineWidth',2);
xlabel('Iteration');
ylabel('Cost');
grid on;

BestSol.Position=pop(13).Position;
figure;
for i = 1:2:numel(BestSol.Position)
    plot (BestSol.Position(1,i) , BestSol.Position(1,i+1),'ro');
    hold on;
    %text (BestSol.Position(1,i) , BestSol.Position(1,i+1), num2str(i/2+0.5),'FontSize',25);
    hold on;
    viscircles ([BestSol.Position(1,i) BestSol.Position(1,i+1)],Rs,'Color', 'k');
end
text (BestSol.Position(1,2) , BestSol.Position(1,3), num2str(1/2+0.5),'FontSize',25);

for i = 1:2:numel(BestSol.Position)
    for j = 1:2:numel(BestSol.Position)
        dist = sqrt((BestSol.Position(1,i)-BestSol.Position(1,j))^2+(BestSol.Position(1,i+1)-BestSol.Position(1,j+1))^2);
        if dist <= Rc && dist ~=0
            plot([BestSol.Position(1,i),BestSol.Position(1,j)],[BestSol.Position(1,i+1),BestSol.Position(1,j+1)],'linewidth',2);
        end
    end
end

grid on;
%}