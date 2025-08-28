%ESTIMATE CONTROLLER PARAMETERS FROM DATA
%load patient data
load("patient_data\pat_2_20_breath_0904_0035.mat")

%EXTRACT DATA
flow = dataSignals.flow; %flow signal
paw = dataSignals.pres; %airway pressure signal
phase= dataSignals.phase;
vol = dataSignals.vol;
time = [0:0.01:length(flow)/100-0.01];

%determine flow as derivative of volume
deriv_vol = diff(vol)/0.01; %forward difference approximation

%plot signals extracted
figure()
plot(flow)
hold on
plot(deriv_vol)
yyaxis right
plot(paw)
hold off

%Process Data
%determine inspiration, determine which points to consider between (start
%and end of a single inspiration)
ind_start = find(diff(phase) == -32);
ind_end = find(diff(phase) == 32);

%extract data between start and end point
%X_array = zeros(20,4);
X_array = zeros(20,4);
X_lin_array = zeros(20,2);

for i = 1:20
    flow_insp = flow(ind_start(i):ind_end(i));
    pres_insp = paw(ind_start(i):ind_end(i));
    time_insp = time(ind_start(i):ind_end(i));
    % 
    % ind_target = find(pres_insp >= 16,1)+6;
    % % 
    % flow_insp = flow_insp(1:ind_target);
    % pres_insp = pres_insp(1:ind_target);
    % time_insp = time_insp(1:ind_target);
    % flow_insp = flow(ind_start(i):ind_start(i)+45);
    % pres_insp = paw(ind_start(i):ind_start(i)+45);
    % time_insp = time(ind_start(i):ind_start(i)+45);
    % 
    ind_shift = -6;
    dot_V = circshift(flow_insp,ind_shift);
    dot_V = dot_V(1:end+ind_shift);
    pres_insp = pres_insp(1:end+ind_shift);

    t = [0:0.01:length(dot_V)/100-0.01]; 
    
    %save in structure
    data_insp.flow = dot_V;
    data_insp.pres = pres_insp;
    data_insp.time = t;
    data_insp.PS = 10;
    data_insp.PEEP = 6;
    data_insp.t = t;
    
    %initial parameters determine from lsqlin
    %create ramp up set point
    int = data_insp.PEEP; %y-intercept of ramp set point signal
    TinspRise = 0.15; %inspiration rise time
    grad = ((data_insp.PS+data_insp.PEEP) - data_insp.PEEP)/TinspRise; %determine gradient of ramp set point sifnal
    
    P_hat = min(data_insp.PS + data_insp.PEEP,grad*data_insp.t + int); %create ramp set point signal
    
    %filter time domain
   % filter =  1/X(1)*exp(-data_insp.t/X(1));
    
    prop_error = P_hat- data_insp.pres';
    % prop_error = (P_hat.*filter) - data_insp.pres';
    int_error = cumsum(prop_error)*0.01;
    
    C = [prop_error',int_error'];
    d = dot_V;
    x_lin = lsqlin(C,d,[],[],[],[],[0,0]);
    X_lin_array(i,:) = x_lin;
    %define ga opt function
    x_hat = [0.0456,4.77e-5,2];
    %x_hat = [0.01,x_hat(1),x_lin(1)/x_lin(2)];
    error_int = ga_control_param(x_hat,data_insp);

    func = @(x_hat)ga_control_param(x_hat,data_insp);
    
    %ga options
    % options = optimoptions('ga', ...
    %   'Display','iter', ...
    %   'FunctionTolerance',1e-6, ...
    %   'PopulationSize', 200, ...
    %   'MaxGenerations', 1000, ...
    %   'UseParallel',true,...
    %   'Display','diagnose');
      
     % 'MaxFunctionEvaluations',50000);
         %'MaxIterations',10000,...

     options = optimoptions('fmincon',...
         'Display','off');

     lb = [1e-9,1e-9,1e-9];
    %lb = [1e-9,1e-9,1e-9];
    %lb = [];
    
    %ub = [1,100,100];
    ub =[];
    %solve genetic algorithm
 
    [X,fval] = fmincon(func,x_hat,[],[],[],[],lb,ub,[],options)
    X_array(i,:) = [X,fval];
    %analyse solution
    %filter time domain
    % filter =  1/X(1)*exp(-data_insp.t/X(1))* 0.01;
   % prop_error = P_hat - data_insp.pres';
    % P_filt =  conv(P_hat, filter, 'same');
    % prop_error = (P_filt.*filter) - data_insp.pres';
   % int_error = cumsum(prop_error)*0.01;
    
    %predicted flow signal
    %dot_V_hat = X(2)*prop_error+(X(2)/X(3))*int_error;
    %dot_V_hat = X(1)*prop_error+(X(1)/X(2))*int_error;
    %measured flow
    
    % figure()
    % plot(dot_V)
    % hold on
    % plot(dot_V_hat)
    % yyaxis right
    % plot(pres_insp)
    % hold on
    % plot(P_hat)
    % plot(P_filt)
    % legend('V_{dot}','V_{dot, hat}','Paw')
    % hold off

    %
    % flow_insp = flow(ind_start(i):ind_end(i));
    % pres_insp = paw(ind_start(i):ind_end(i));
    % ind_shift = -6;
    % dot_V = circshift(flow_insp,ind_shift);
    % dot_V = dot_V(1:end+ind_shift);
    % pres_insp = pres_insp(1:end+ind_shift);
    % t = [0:0.01:length(dot_V)/100-0.01]; 
    % P_hat = min(data_insp.PS + data_insp.PEEP,grad*t + int); %create ramp set point signal
    % prop_error = P_hat - pres_insp';
    % % P_filt =  conv(P_hat, filter, 'same');
    % % prop_error = (P_filt.*filter) - data_insp.pres';
    % int_error = cumsum(prop_error)*0.01;
    % 
    % %predicted flow signal
    % %dot_V_hat = X(2)*prop_error+(X(2)/X(3))*int_error;
    % dot_V_hat = X(1)*prop_error+(X(1)/X(2))*int_error;
    % 
   
    Kp = X(1);
    Ti = X(2);
    T1 = X(3);
    % T2 = 1;
    % T1 = T2*10;
    % sys = Kp*(T1*s+1)/(T2*s+1);
 
    controller_1 = Kp*(1+1/(Ti*s));
    controller_2 = Kp/Ti*(Ti*s+1)/(s+X(3));
    sys = controller_2;%*(1/(1+tau*s));
    
    Y = lsim(sys,prop_error,t);

    figure()
    plot(t,Y)
    hold on
    plot(t,dot_V)
    yyaxis right
    plot(t, prop_error)
    legend("dot_v_hat","dot_v", "error")
    hold off



    % figure()
    % plot(dot_V)
    % hold on
    % plot(dot_V_hat)
    % yyaxis right
    % plot(pres_insp)
    % hold on
    % % plot(prop_error)
    % % plot(int_error)
    % % plot(X(1)/X(2)*int_error)
    % plot(P_hat)
    % plot(P_filt)
    % legend('V_{dot}','V_{dot, hat}','Paw')
    % hold off
    
end

function error_opt = ga_control_param(x_hat,data_insp)
    %tau = x_hat(1)
    %K_p = x_hat(2)
    %T_i = x_hat(3)

    %create ramp up set point
    int = data_insp.PEEP; %y-intercept of ramp set point signal
    TinspRise = 0.15; %inspiration rise time
    grad = ((data_insp.PS+data_insp.PEEP) - data_insp.PEEP)/TinspRise; %determine gradient of ramp set point sifnal
    
    P_hat = min(data_insp.PS + data_insp.PEEP,grad*data_insp.t + int); %create ramp set point signal
    
    %filter time domain
    %filter =  1/x_hat(1)*exp(-data_insp.t/x_hat(1))* 0.01;
    prop_error = P_hat - data_insp.pres';
    s = tf('s');
    controller_2 = x_hat(1)/x_hat(2)*(x_hat(2)*s+1)/(s+x_hat(3));
    sys = controller_2;%*(1/(1+tau*s));
    dot_V_hat = lsim(sys,prop_error,data_insp.t);

    %P_filt =  conv(P_hat, filter, 'same') ;
   % prop_error = (P_filt.*filter) - data_insp.pres';
    % int_error = cumsum(prop_error)*0.01;
     
    %predicted flow signal
    %dot_V_hat = x_hat(2)*prop_error+(x_hat(2)/x_hat(3))*int_error;
    % dot_V_hat = x_hat(1)*prop_error+(x_hat(1)/x_hat(2))*int_error;
    %measured flow
    %ind_shift = -6;
    %dot_V = circshift(data_insp.flow,ind_shift);
    dot_V = data_insp.flow;

    % figure()
    % plot(dot_V)
    % hold on
    % plot(dot_V_hat)
    
    %error
    error_opt = sum((dot_V-dot_V_hat).^2);
end
