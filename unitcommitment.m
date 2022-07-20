% Unit Commitment by Yalmip
clear
%% Part 1 Simple Example
% Data Setting
N_unit = 3;
Horizon = 48;

P_max = [100;50;25];
P_min = [20;40;1];

% Cost of Unit
Q = diag([.04 .01 .02]); %Quatric
C = [10 20 20];

P_forecast = 100 + 50*sin((1:Horizon)*2*pi/24); % 2 Period

% Variables
Onoff = binvar(N_unit,Horizon,'full');
P = sdpvar(N_unit,Horizon,'full');

% Basic Constraints
Cons = [];

% unit output limit
for k=1:Horizon
    Cons = [Cons, Onoff(:,k).*P_min <= P(:,k) <= Onoff(:,k).*P_max];
end
for k=1:Horizon
    Cons = [Cons, sum(P(:,k)) >= P_forecast(k)];
end

% Basic Objective
Obj = 0;
for k=1:Horizon
    Obj = Obj + P(:,k)' * Q * P(:,k) + C * P(:,k);
end

% Optimize Setting
Ops = sdpsettings('verbose', 1, 'debug', 1);
optimize(Cons, Obj, Ops);
clf % clear figure
stairs(value(P)');
legend('Unit 1','Unit 2','Unit 3')

%% Part 2 min up- & down- time
minup   = [6;30;1];
mindown = [3;6;3];

% turn on
for k=2:Horizon
    for unit=1:N_unit
        % indicator will be 1 only when turn on
        indicator = Onoff(unit,k) - Onoff(unit,k-1); 
        range = k:min(Horizon,k+minup(unit)-1);
        % add constraints
        Cons = [Cons, Onoff(unit,range) >= indicator];
    end
end

% turn down
for k=2:Horizon
    for unit=1:N_unit
        % indicator will be 1 only when turn down
        indicator = Onoff(unit,k-1) - Onoff(unit,k); 
        range = k:min(Horizon,k+mindown(unit)-1);
        % add constraints
        Cons = [Cons, Onoff(unit,range) <= 1-indicator];
    end
end

% optimize
Ops = sdpsettings('verbose', 2, 'debug', 1);
optimize(Cons, Obj, Ops);
clf % clear figure
stairs(value(P)');
legend('Unit 1','Unit 2','Unit 3')

%% Part 3 Quantized power-level
Unit3_level = [0 1 6 10 12 20];
for k=1:Horizon
    Cons = [Cons, ismember(P(3,k),Unit3_level)];
end
optimize(Cons, Obj, Ops);
clf % clear figure
stairs(value(P)');
legend('Unit 1','Unit 2','Unit 3')

%% Part 4 Efficient simulation
% REMARKS: 
% real output of origin code is different from tutorial since the Change Penalty
% is too high and the power balance equation is generation > demand which
% made the optimal output of units stay at high generate level
clear
% Basic Settings
N_unit = 3;
Horizon = 48;

P_max = [100;50;25];
P_min = [20;40;1];

% Cost of Unit
Q = diag([.04 .01 .02]); %Quatric
C = [10 20 20];

P_forecast = 100 + 50*sin((1:Horizon)*2*pi/24); % 2 Period

minup   = [6;30;1];
mindown = [3;6;3];
% Variables
Onoff = binvar(N_unit,Horizon,'full');
P = sdpvar(N_unit,Horizon,'full');

N_history = max([minup;mindown]);
P_forecast = sdpvar(1,Horizon);
HistoryOnoff = binvar(N_unit,N_history);
DemandSlack = sdpvar(1,Horizon);
Pre_P = sdpvar(N_unit,1);
DemandPenalty = 1000;
ChangePenalty = 10;
% optimize setting
Cons = [];
Obj = 0;
% Optimize Basic Settings
for k=1:Horizon
    % unit output constraints
    Cons = [Cons, Onoff(:,k).*P_min <= P(:,k) <= Onoff(:,k).*P_max];
    % power balance
    Cons = [Cons, sum(P(:,k)) + DemandSlack(k) >= P_forecast(:,k)];
    % demand slack constraints
    Cons = [Cons, DemandSlack(:,k) >= 0]
    % Objective
    Obj = Obj + P(:,k)' * Q * P(:,k) + C * P(:,k);
    Obj = Obj + DemandSlack(k)*DemandPenalty;
end
% add on-off constraints
Cons = [Cons, cons_equtive_ON([HistoryOnoff Onoff], minup)];
Cons = [Cons, cons_equtive_ON(1-[HistoryOnoff Onoff], mindown)];
% add control-change penalty
for k=2:Horizon
    Obj = Obj + ChangePenalty * norm(P(:,k)-P(:,k-1),1);
end
Obj = Obj + ChangePenalty * norm(P(:,1)-Pre_P,1);

% create OPTIMIZER object
Parameters = {HistoryOnoff, P_forecast, Pre_P};
Outputs = {P, Onoff};
ops = sdpsettings('verbose',2,'debug',1);
% OPTIMIZER!!!
Controller = optimizer(Cons,Obj,ops,Parameters,Outputs);
% optimizer creates an object with a pre-compiled low-level numerical format 
% which can be used to efficiently solve a series of similar problems 
% (reduces YALMP analysis and compilation overhead)

% create history data
OldOnoff = repmat([1;1;0],1,N_history);
OldP = repmat([100;40;0],1,N_history);

% *********** SIMULATION **************
for k=1:500
    % forecast baseline
    forecast = 100;
    % Daily fluctation
    forecast = forecast + 50*sin((k:k+Horizon-1)*2*pi/24); % window = 48
    % Add other effect
    forecast = forecast + 20*sin((k:k+Horizon-1)*2*pi/24/7);
    % Add disturbance
    forecast = forecast + randn(1,Horizon)*5;
    
    [solution,problem] = Controller{OldOnoff,forecast,OldP(:,end)};
    optimalP = solution{1};
    optimalOnoff = solution{2};
    
    hold off
    stairs([-N_history+1:0 1:Horizon],[OldP optimalP]');
    hold on
    stairs(1:Horizon,forecast,'k+'); % mark load demand with black +
    axis([-N_history Horizon 0 170]);
    drawnow
    pause(0.05);
    % shift history data
    OldP = [OldP(:,2:end) optimalP(:,1)];
    OldOnoff = [OldOnoff(:,2:end) optimalOnoff(:,1)];
end


% Add turn-on -off constraints function
function Con = cons_equtive_ON(x,minup)
    if min(size(x))==1 % only 1 time point
        x = x(:)'; 
    end
    if size(x,1) ~= size(minup,1)
        error('MINUP should have as many rows as X');
    end
    time_point = size(x,2);
    Con = [];
    for k=2:time_point
        for unit=1:size(x,1)
            % indicator will be 1 only when turn on
            indicator = x(unit,k) - x(unit,k-1); 
            range = k:min(time_point,k+minup(unit)-1);
            % add Constraints
            affected = x(unit,range);
            if strcmp(class(affected),'binvar')
                Con = [Con,affected >= indicator];
            end
        end
    end
end


