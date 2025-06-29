%
% Copyright (c) 2015, Mostapha Kalami Heris & Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "LICENSE" file for license terms.
%
% Project Code: YPEA124
% Project Title: Implementation of MOEA/D
% Muti-Objective Evolutionary Algorithm based on Decomposition
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Cite as:
% Mostapha Kalami Heris, MOEA/D in MATLAB (URL: https://yarpiz.com/95/ypea124-moead), Yarpiz, 2015.
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

clc;
clear;
close all;

%% Problem Definition
VarMax = 100;
VarMin = 0;
Covered_Area = zeros(100,100);
Obstacle_Area = ones(100,100);
N = 60;
rc = 30; 
rs = ones(1,N)*10;
stat(1,:)=rs;
stat(2,1)=rc;
sink=[50 50];

% Number of Objective Functions
nObj = 2;


%% MOEA/D Settings

MaxIt = 1000;  % Maximum Number of Iterations

nPop = 100;    % Population Size (Number of Sub-Problems)

nArchive = 50;

T = max(ceil(0.15*nPop), 2);    % Number of Neighbors
T = min(max(T, 2), 15);

crossover_params.gamma = 0.5;
crossover_params.VarMin = VarMin;
crossover_params.VarMax = VarMax;

%% Initialization

% Create Sub-problems
sp = CreateSubProblems(nObj, nPop, T);

% Empty Individual
empty_individual.Position = [];
empty_individual.Cost = [];
empty_individual.g = [];
empty_individual.IsDominated = [];

% Initialize Goal Point
z = inf(nObj, 1);
%z = zeros(nObj, 1);

% Create Initial Population
pop = repmat(empty_individual, nPop, 1);
for i = 1:nPop
    pop(i).Position = unifrnd(0,100,[N 2]);
    pop(i).Position(1,:) = sink;
    pop(i).Cost = CostFunction(pop(i).Position,stat,Obstacle_Area,Covered_Area);
    z = min(z, pop(i).Cost);
end

for i = 1:nPop
    pop(i).g = DecomposedCost(pop(i), z, sp(i).lambda);
end

% Determine Population Domination Status
pop = DetermineDomination(pop);

% Initialize Estimated Pareto Front
EP = pop(~[pop.IsDominated]);

%% Main Loop

for it = 1:MaxIt
    for i = 1:nPop
        
        % Reproduction (Crossover)
        K = randsample(T, 2);
        
        j1 = sp(i).Neighbors(K(1));
        p1 = pop(j1);
        
        j2 = sp(i).Neighbors(K(2));
        p2 = pop(j2);
        
        y = empty_individual;
        y.Position = Crossover(p1.Position, p2.Position, crossover_params);
        y.Position(1,:) = sink;
        
        y.Cost = CostFunction(y.Position,stat,Obstacle_Area,Covered_Area);
        
        z = min(z, y.Cost);
        
        %for j = sp(i).Neighbors
        for j = sp(i).Neighbors(1:1)
            y.g = DecomposedCost(y, z, sp(j).lambda);
            if y.g < pop(j).g
                pop(j) = y;
            end
        end
        
    end
    
    % Determine Population Domination Status
	pop = DetermineDomination(pop);
    
    ndpop = pop(~[pop.IsDominated]);
    
    EP = [EP
        ndpop]; %#ok
    
    EP = DetermineDomination(EP);
    EP = GetParetoFront(EP);
    %EP = EP(~[EP.IsDominated]);
    
    if numel(EP)>nArchive
        Extra = numel(EP)-nArchive;
        ToBeDeleted = randsample(numel(EP), Extra);
        EP(ToBeDeleted) = [];
    end
    
    % Plot EP
    figure(1);
    PlotCosts(EP);
    pause(0.01);

    
    % Display Iteration Information
    disp(['Iteration ' num2str(it) ': Number of Pareto Solutions = ' num2str(numel(EP))]);
    
end

%% Reults

disp(' ');

EPC = [EP.Cost];
for j = 1:nObj
    
    disp(['Objective #' num2str(j) ':']);
    disp(['      Min = ' num2str(min(EPC(j, :)))]);
    disp(['      Max = ' num2str(max(EPC(j, :)))]);
    disp(['    Range = ' num2str(max(EPC(j, :))-min(EPC(j, :)))]);
    disp(['    St.D. = ' num2str(std(EPC(j, :)))]);
    disp(['     Mean = ' num2str(mean(EPC(j, :)))]);
    disp(' ');
    
end


