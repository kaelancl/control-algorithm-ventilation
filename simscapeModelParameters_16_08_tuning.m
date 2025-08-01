%Parameters for Simscape Model
clear all
tic
simTime = 2002;
t = [0:0.01:(simTime-0.01)];

%load('PmusHatSmoothNorm_30_10.mat')

%PmusHatSmoothNorm = PmusHatSmoothNorm(1:end-1);
%plot(PmusHatSmoothNorm)
%hold off

%VENTILATOR SETTINGS
PEEP = 5; %cmH2O, Posistive End Expiratory Pressure
Paw_Ref = 11; %Pressure target, cmH_2O

%Reference support target changes after 60s
Paw_Ref_in = [t', ones(length(t),1)*Paw_Ref];
Paw_Ref_in(7000:end,2) = Paw_Ref+2;
PS = Paw_Ref - PEEP;

maxAllowableChange = 5; %+- max allowed change in pressure support (cmH2O)
Paw_Threshold = -1; %cmH2O, trigger variable
Q_Threshold_Percentage = 0.36; %percent of peak flow, cycling off variable
numOfBreaths = 20; %number of tidal breaths until support adjust

timeInspVent = 2;
timeExpVent = 2;

cutInValue = 13;
alpha = -1/(cutInValue-15);
zeta = 0.3258;
beta = 3.2576;

%CONTROLLER
%PID
Kp = 3; %Proportional constant, PID cntrol
Ti = 0.1; %Integral time constant, PID control
D = 0; %Derivative constant, PID control

P = 1;
I = 1;
D = 0;


%Adjustment Support Level Filter
timeResponse = 25; %time response of system, amount of time for controller
%adjust support level

%Muscle Pressure Adjustment Curve Parameters
gradPmus = 2.5; %gradient of muscle pressure support level adjustment curve
yIntPmus = -20; %y intercept
a1 = 0.1; %gain for driving pressure driven support level adjustment
%Driving Pressure Adjustment
gradDrivePres = 2.5;
yIntDrivePres = -32.5;
a2 = 0.1; %gain for muscle pressure driven support level adjustment

%INITIAL PATIENT SETTINGS
Resistance = 8.25; %cmH2O/L/s
Compliance = 1/13.8; %L/cmH2O
%initialMaxEffort = 10.5; %-cmH2O, initial maximum effort of patient
%maximum effort changes as PSV level changes
timeInsp = 1.4; %patient inspiration time, seconds
timeExp = 2.6; %patient expiration, seconds

varPE = 10; % %change in variation effort (relative to max effort)

%OCCLUSION SETTINGS
%end expiratory occlusion
lengthEEOcc = 5; %Length of end-expiratory occlusion (seconds)
breathsTillEE = 3; %number of breaths that occur before EE occlusion
numberOfRecordingsEE = 3;

%end inspiratory occlusion
numberOfRecordingsEI = 3; %number of recordings till support level adjusted
lengthEIOcc = 1.5; %Length of end-inspiratory occlusion (seconds)
flowOccThreshold = 0.05; %percentage of peak flow when occlusion triggered
breathsTillEI = 3; %number of breaths that occur before EI occlusion


%state (and inverse state signal)
tCounter = 0;
state = zeros(1,length(tCounter));
invState = zeros(1,length(tCounter));
mode = zeros(1,length(tCounter));
breathCounter = 0;
x = 0;
Pmus = zeros(1,length(tCounter));
flowIndexEnable = zeros(1,length(tCounter));

for idx=1:length(t)
    if tCounter < timeInspVent
        state(idx) = 1;
        invState(idx) = 0;
        mode(idx) = 1;
        flowIndexEnable(idx) = 1;
    end
    %expiration
    if breathCounter == 100000
        if tCounter >= timeInspVent && tCounter < (timeInspVent + lengthEIOcc)
            state(idx) = 0;
            invState(idx) = 0;
            mode(idx) = 2;
            
              if tCounter < (timeInspVent + 0.1)
                flowIndexEnable(idx) = 1;
              else
                  flowIndexEnable(idx) = 0;
              end
        end
        if tCounter > (timeInspVent + lengthEIOcc) && tCounter < (timeInspVent + timeExpVent + lengthEIOcc)
            state(idx) = 0;
            invState(idx) = 1;
            mode(idx) = 0;
        end
        if tCounter >= (timeInspVent + timeExpVent + lengthEIOcc)
            tCounter = 0;
            breathCounter = 0;
        end
    else
        if tCounter >= timeInspVent && tCounter < (timeInspVent + timeExpVent)
            state(idx) = 0;
            invState(idx) = 1;
            mode(idx) = 0;
            if tCounter < (timeInspVent + 0.1)
                flowIndexEnable(idx) = 1;
            else 
                flowIndexEnable(idx) = 0;
            end
        end
         if tCounter >= (timeInspVent + timeExpVent)
            tCounter = 0;
            breathCounter = breathCounter + 1;
        end
    end
  
    tCounter = tCounter + 0.01;
end




%plot(t(1:2000),-Pmus(1:2000))
effortMode = 0;
if effortMode == 1
    initialMaxEffort = zeros(1, length(t));
    x = 0;
    for idx = 1:length(t)

        if x < 7
            initialMaxEffort(idx) = 2;
        elseif x < 14
            initialMaxEffort(idx) = 2.25;
        elseif x < 21
            initialMaxEffort(idx) = 2.5;
        elseif x < 28
            initialMaxEffort(idx) = 2.75;
        elseif x < 35
            initialMaxEffort(idx) = 3;
        elseif x < 42
            initialMaxEffort(idx) = 3.25;
        elseif x < 49
            initialMaxEffort(idx) = 3.5;
        elseif x < 56
            initialMaxEffort(idx) = 3.75;
        elseif x < 63
            initialMaxEffort(idx) = 4;
        elseif x < 70
            initialMaxEffort(idx) = 4.25;
        elseif x < 77
            initialMaxEffort(idx) = 4.5;
        elseif x < 84
            initialMaxEffort(idx) = 4.75;
        elseif x < 91
            initialMaxEffort(idx) = 5;
        elseif x < 98
            initialMaxEffort(idx) = 5.25;
        elseif x < 105
            initialMaxEffort(idx) = 5.5;
        elseif x < 112
            initialMaxEffort(idx) = 5.75;
        elseif x < 119
            initialMaxEffort(idx) = 6;
        elseif x < 126
            initialMaxEffort(idx) = 6.25;
        elseif x < 133
            initialMaxEffort(idx) = 6.5;
        else
        end
        % elseif x < 1100
        %     initialMaxEffort(idx) = 5.5;
        % elseif x < 1200
        %     initialMaxEffort(idx) = 6;
        % elseif x < 1300
        %     initialMaxEffort(idx) = 6.5;
        % elseif x < 1400
        %     initialMaxEffort(idx) = 7;
        % elseif x < 1500
        %     initialMaxEffort(idx) = 7.5;
        % elseif x < 1600
        %     initialMaxEffort(idx) = 8;
        % elseif x < 1700
        %     initialMaxEffort(idx) = 8.5;
        % elseif x < 1800
        %     initialMaxEffort(idx) = 9;
        % elseif x < 1900
        %     initialMaxEffort(idx) = 9.5;
        % elseif x < 2000
        %     initialMaxEffort(idx) = 10;
        % elseif x < 2100
        %     initialMaxEffort(idx) = 10.5;
        % elseif x < 2200
        %     initialMaxEffort(idx) = 11;
        % elseif x < 2300
        %     initialMaxEffort(idx) = 11.5;
        % elseif x < 2400
        %     initialMaxEffort(idx) = 12;
        % elseif x < 2500
        %     initialMaxEffort(idx) = 12.5;
        % elseif x < 2600
        %     initialMaxEffort(idx) = 13;
        % elseif x < 2700
        %     initialMaxEffort(idx) = 13.5;
        % elseif x < 2800
        %     initialMaxEffort(idx) = 14;
        % elseif x < 2900
        %     initialMaxEffort(idx) = 14.5;
        % elseif x < 3000
        %     initialMaxEffort(idx) = 15;
        % end
        x = x + 0.01;
    end
    effortIn = [t' initialMaxEffort'];
else
    initialMaxEffort = 10.35+2;
end


%% Pmus Model
load('general_model.mat')

%numBreaths = simTime/((length(pmus_model)-1)/100);
%lengthBreath = length(pmus_model)-1;
%startPos = 1;
Pmus = [];
for idx=1:simTime/7
    Pmus = [Pmus pmus_model];
    %startPos = startPos + lengthBreath;
    %Pmus(startPos+1:startPos+231) = 0;
    %startPos = startPos + 231;
end


stateIn = [t',state'];
invStateIn = [t',invState'];
Pmus_in = [t', Pmus(1:end)'];
modeIn = [t', mode'];

%% Patient Response Variables

m = -1.22;
ratioPV = ones(1, length(t))*m;
ratioPV_in = [t', ratioPV'];
c = 23.56;


% C_prop = 0; %propofol concentration [ug*mL^-1]
% S_0 = 0; %receptor sensitivity in hyperoxia [L*min^-1*(nM*L^-1)^-1]
% A = 17.8; %Area constant
% P_o2 = 100; %partial pressure of oxygen
% P_0 = 30; %P_o2 max sensitivity beffore receptor failure
% T_c = 40; %central receptors threshold [nM L-1]
% T_p = 34.6; %peripheral receptor threshold [nM L-1]
% F_R = 12; %respiratory rate [min^-1]
% V_D_anat = 0.12; %anatomical dead space [L]
% V_D_alv_ratio = 0.3; %alveolar dead space to tidal volume ratio
% dot_V_TH = 0.2689; %flow threshold [L*s^-1]
% D_w = 0; %wakefulness drive (0 = asleep)
% Na_conc = 0.139;
% K_conc = 0.0043;
% Ca_conc = 0.0012;
% Mg_conc = 0.00085;
% Alb = 4.5;
% Pi = 0.00012;
% HCO3_conc = 0.024; %nM*L-1
% V_co2 = 0.2; %L * min^-1
% T_i = 0.978; %Inspiration time [s]
% Cl_conc = 0.105; %Chlorine concentration [M*L^-1]
% K_w = 2.39e-14; %ion product of water 
% K_c = 2.45e-11; %equilibrium constant
% K_2 = 2.19e-7; %phosphoric acid dissociation constant
% K_3 = 1.16e-10; %carbonate dissocation constant
% K_h = 1.77e-7; %histidine dissociation constant
% 
% PEEP = 6; %cmH2O, Posistive End Expiratory Pressure
% Resistance = 11.2; %cmH2O/L/s
%Compliance = 0.0617; %L/cmH2O
% figure()
% plot(t, mode')
% 
% figure()
% plot(t, state')
% 
% figure()
% plot(t, invState')
% 
% figure()
% plot(t, Pmus')
% plot(t, state)
% hold on
% plot(t, invState)
% plot(t, Pmus)
% plot(t, mode)

% plot(modeIn, 'Color', 'r')    
% hold on
% plot(flowIndexEnableIn, 'Color', 'g')
% hold off

%SIMULATE MODEL
%simOut = sim("simplePatientModelTuning_30_08_23.slx");
% simOut = sim("simScapeModel_tuning_08_23.slx")
% %DETERMINE IAE
% %take average of last set point values from 380 seconds onwards
% % ssIndStart = find(simOut.adjustedSetPoint.time > 400,1);
% % ssValue = mean(simOut.adjustedSetPoint.data(ssIndStart:end));
% % 
% % changeFlagTime = simOut.changeFlag.time(find(simOut.changeFlag.data ~= 0));
% % changeFlagTime = round(changeFlagTime*100)/100;
% % changeFlagTime = unique(changeFlagTime);
% % for idx=2:length(changeFlagTime)
% %     if abs(diff([changeFlagTime(idx) changeFlagTime(idx-1)])) < 0.1
% %         changeFlagTime(idx) = 0;
% %     end
% % end
% % changeFlagTimeProcessed = changeFlagTime(find(changeFlagTime ~= 0));
% % IAE = 0;
% % IAE_array = [];
% % 
% % %DOESN"T QUITE WORK, I THINK changeFlagTime is a bit funky
% % % for idx = 1:length(changeFlagTimeProcessed)
% % %     IAE = IAE + abs(diff([ssValue simOut.adjustedSetPoint.data(find(simOut.adjustedSetPoint.time > changeFlagTimeProcessed(idx),1))]));
% % %     IAE_array = [IAE_array; IAE];
% % % end
% % 
% % for idx = 1:length(changeFlagTimeProcessed)
% %     IAE = IAE + abs(diff([ssValue simOut.adjustedSetPoint.data(find(simOut.adjustedSetPoint.time > changeFlagTimeProcessed(idx),1))]));
% %     IAE_array = [IAE_array; IAE];
% % end
% % 
% plot(simOut.adjustedSetPoint.time,simOut.adjustedSetPoint.data)
% hold on
% % plot(changeFlagTimeProcessed, IAE_array)
% % toc
% 
% %SETTLING TIME
% PSVout= simOut.adjustedSetPoint.data;
% TimeOut = simOut.adjustedSetPoint.time;
% ssValue = PSVout(end);
% upper = yline(ssValue + (0.02*ssValue))
% lower = yline(ssValue - (0.02*ssValue))
% diffFromSS = abs(PSVout - ssValue);
% ssBand = ssValue*0.02;
% settlingTime = TimeOut(find(diffFromSS > ssBand,1, 'last')) 
% settlingPoint = PSVout(find(TimeOut == settlingTime))
% plot(settlingTime, settlingPoint, '-o')
% %settlingTimeIdx = find(abs(diff([reversePSV(1) reversePSV])) > 0.05*reversePSV(1))


% selectPmus = 1;
% %Base Pmus signal
% breathCounter = 0;
% cycleCounter = 1;
% a = [1.2 1.5 2 3 6];
% b = [6 3 2 1.5 1.2];
% c = [0.05:0.05:0.75];
% d = [0.1:0.1:1.5];
% e = [2:0.5:6];
% 
% switch selectPmus
%     case 1 %half sine wave    
%         for idx = 1:length(t)
%             if x <= timeInsp %patient effort during inspiration
%                 Pmus(idx) = sin(x/timeInsp*pi);
%             end
% 
%             if breathCounter == 3
%                 if x > timeInsp && x <= (timeInsp + timeExp + lengthEIOcc) %patient effort during expiration
%                     Pmus(idx) = 0;
%                 end
% 
%                 if x > (timeInsp + timeExp + lengthEIOcc)%reset timer after 2 seconds has passed
%                     x = 0;
%                     breathCounter = 0;
%                 end
%             else
%                 if x > timeInsp && x <= timeInsp + timeExp %patient effort during expiration
%                     Pmus(idx) = 0;
%                 end
% 
%                 if x > timeInsp + timeExp %reset timer after 2 seconds has passed
%                     x = 0;
%                     breathCounter = breathCounter + 1;
%                 end
%             end
%             x = x + 0.01; %increment timer
%         end
% 
%     case 2 %flipped halfwave
% 
%         Tip = 0.5;
%         Te = 1.5;
%         Tep = 2;
%         Tnext = 4;
% 
%         for idx = 1:length(t)
%             if x <Tip
%                 Pmus(idx) = -sin(((x)*pi)/(2*(-Tip)));
%             end
%             if  x >= Tip && x < Te
%                 Pmus(idx) = 1;
%             end
% 
%              if x >= Te && x < Tep
%                 Pmus(idx) = sin(((x-Te)*pi)/(2*(Tep-Te))+pi)+1 ;
%              end         
%             if breathCounter == 3
%                 if x > Tep && x <= (Tnext + lengthEIOcc) %patient effort during expiration
%                     Pmus(idx) = 0;
%                 end
% 
%                 if x > (Tnext + lengthEIOcc)%reset timer after 2 seconds has passed
%                     x = 0;
%                     breathCounter = 0;
%                 end
%             else
%                  if x >=Tep && x < Tnext
%                     Pmus(idx) = 0;
%                  end
%                 if x >= Tnext
%                     x = 0;
%                     breathCounter = breathCounter + 1;
%                 end
%             end
%             x = x + 0.01;
%         end
% 
%     case 3 %square wave    
%         for idx = 1:length(t)
%             if x <= timeInsp %patient effort during inspiration
%                 Pmus(idx) = 1;
%             end
% 
%             if breathCounter == 3
%                 if x > timeInsp && x <= (timeInsp + timeExp + lengthEIOcc) %patient effort during expiration
%                     Pmus(idx) = 0;
%                 end
% 
%                 if x > (timeInsp + timeExp + lengthEIOcc)%reset timer after 2 seconds has passed
%                     x = 0;
%                     breathCounter = 0;
%                 end
%             else
%                 if x > timeInsp && x <= timeInsp + timeExp %patient effort during expiration
%                     Pmus(idx) = 0;
%                 end
% 
%                 if x > timeInsp + timeExp %reset timer after 2 seconds has passed
%                     x = 0;
%                     breathCounter = breathCounter + 1;
%                 end
%             end
%             x = x + 0.01; %increment timer
%         end
% 
%     case 4 %saw tooth
%         for idx = 1:length(t)
%             if x <= timeInsp %patient effort during inspiration
%                 if x<= timeInsp/2
%                      Pmus(idx) = 4/3*x;
%                 elseif x > timeInsp/2
%                     Pmus(idx) = -4/3*x+2;
%                 end
%             end
% 
%             if breathCounter == 3
%                 if x > timeInsp && x <= (timeInsp + timeExp + lengthEIOcc) %patient effort during expiration
%                     Pmus(idx) = 0;
%                 end
% 
%                 if x > (timeInsp + timeExp + lengthEIOcc)%reset timer after 2 seconds has passed
%                     x = 0;
%                     breathCounter = 0;
%                 end
%             else
%                 if x > timeInsp && x <= timeInsp + timeExp %patient effort during expiration
%                     Pmus(idx) = 0;
%                 end
% 
%                 if x > timeInsp + timeExp %reset timer after 2 seconds has passed
%                     x = 0;
%                     breathCounter = breathCounter + 1;
%                 end
%             end
%             x = x + 0.01; %increment timer
%         end
% 
%     case 5 %sawtooth shifting peak
%         for idx = 1:length(t)
%             if x <= timeInsp %patient effort during inspiration
%                Pmus(idx) = min(-a(cycleCounter)/timeInsp*x+a(cycleCounter),b(cycleCounter)*x/timeInsp);
% 
%             end
% 
%             if breathCounter == 3
%                 if x > timeInsp && x <= (timeInsp + timeExp + lengthEIOcc) %patient effort during expiration
%                     Pmus(idx) = 0;
%                 end
% 
%                 if x > (timeInsp + timeExp + lengthEIOcc)%reset timer after 2 seconds has passed
%                     x = 0;
%                     breathCounter = 0;
%                     cycleCounter = cycleCounter + 1;
%                     if cycleCounter > 5
%                         cycleCounter = 5;
%                     end
%                 end
%             else
%                 if x > timeInsp && x <= timeInsp + timeExp %patient effort during expiration
%                     Pmus(idx) = 0;
%                 end
% 
%                 if x > timeInsp + timeExp %reset timer after 2 seconds has passed
%                     x = 0;
%                     breathCounter = breathCounter + 1;
%                 end
%             end
%             x = x + 0.01; %increment timer
%         end
% 
%     case 6 %short triangle
%         for idx = 1:length(t)
%             if x <= timeInsp %patient effort during inspiration
% 
%                Pmus(idx) = min(1/0.125*x-5,-x/0.125+7);
%                Pmus(idx) = max(Pmus(idx),0);
% 
%             end
% 
%             if breathCounter == 3
%                 if x > timeInsp && x <= (timeInsp + timeExp + lengthEIOcc) %patient effort during expiration
%                     Pmus(idx) = 0;
%                 end
% 
%                 if x > (timeInsp + timeExp + lengthEIOcc)%reset timer after 2 seconds has passed
%                     x = 0;
%                     breathCounter = 0;
%                     cycleCounter = cycleCounter + 1;
%                     if cycleCounter > 5
%                         cycleCounter = 5;
%                     end
%                 end
%             else
%                 if x > timeInsp && x <= timeInsp + timeExp %patient effort during expiration
%                     Pmus(idx) = 0;
%                 end
% 
%                 if x > timeInsp + timeExp %reset timer after 2 seconds has passed
%                     x = 0;
%                     breathCounter = breathCounter + 1;
%                 end
%             end
%             x = x + 0.01; %increment timer
%         end
%     case 7 %triangle varying lengths
%         for idx = 1:length(t)
%             if x <= timeInsp %patient effort during inspiration
% 
%                Pmus(idx) = min(1/c(cycleCounter)*x,-x/c(cycleCounter)+2);
%                Pmus(idx) = max(Pmus(idx),0);
% 
%             end
% 
%             if breathCounter == 3
%                 if x > timeInsp && x <= (timeInsp + timeExp + lengthEIOcc) %patient effort during expiration
%                     Pmus(idx) = 0;
%                 end
% 
%                 if x > (timeInsp + timeExp + lengthEIOcc)%reset timer after 2 seconds has passed
%                     x = 0;
%                     breathCounter = 0;
%                     cycleCounter = cycleCounter + 1;
%                     if cycleCounter > 15
%                         cycleCounter = 15;
%                     end
%                 end
%             else
%                 if x > timeInsp && x <= timeInsp + timeExp %patient effort during expiration
%                     Pmus(idx) = 0;
%                 end
% 
%                 if x > timeInsp + timeExp %reset timer after 2 seconds has passed
%                     x = 0;
%                     breathCounter = breathCounter + 1;
%                 end
%             end
%             x = x + 0.01; %increment timer
%         end
%     case 8 %sinusoidal wave varying lengths
%         for idx = 1:length(t)
%             if x <= timeInsp %patient effort during inspiration
%                if x < d(cycleCounter)
%                     Pmus(idx) = sin(pi/d(cycleCounter)*x);
%                else 
%                    Pmus(idx) = 0;
%                end
%             end
% 
%             if breathCounter == 3
%                 if x > timeInsp && x <= (timeInsp + timeExp + lengthEIOcc) %patient effort during expiration
%                     Pmus(idx) = 0;
%                 end
% 
%                 if x > (timeInsp + timeExp + lengthEIOcc)%reset timer after 2 seconds has passed
%                     x = 0;
%                     breathCounter = 0;
%                     cycleCounter = cycleCounter + 1;
%                     if cycleCounter > 15
%                         cycleCounter = 15;
%                     end
%                 end
%             else
%                 if x > timeInsp && x <= timeInsp + timeExp %patient effort during expiration
%                     Pmus(idx) = 0;
%                 end
% 
%                 if x > timeInsp + timeExp %reset timer after 2 seconds has passed
%                     x = 0;
%                     breathCounter = breathCounter + 1;
%                 end
%             end
%             x = x + 0.01; %increment timer
%         end
% 
%     case 9 %triangle wave varying lengths, capped at 1
%         for idx = 1:length(t)
%             if x <= timeInsp %patient effort during inspiration
%               Pmus(idx) = min([1/0.25*x,-x/0.25+e(cycleCounter),1]);
%                Pmus(idx) = max(Pmus(idx),0);
%             end
% 
%             if breathCounter == 3
%                 if x > timeInsp && x <= (timeInsp + timeExp + lengthEIOcc) %patient effort during expiration
%                     Pmus(idx) = 0;
%                 end
% 
%                 if x > (timeInsp + timeExp + lengthEIOcc)%reset timer after 2 seconds has passed
%                     x = 0;
%                     breathCounter = 0;
%                     cycleCounter = cycleCounter + 1;
%                     if cycleCounter > 9
%                         cycleCounter = 9;
%                     end
%                 end
%             else
%                 if x > timeInsp && x <= timeInsp + timeExp %patient effort during expiration
%                     Pmus(idx) = 0;
%                 end
% 
%                 if x > timeInsp + timeExp %reset timer after 2 seconds has passed
%                     x = 0;
%                     breathCounter = breathCounter + 1;
%                 end
%             end
%             x = x + 0.01; %increment timer
%         end
%     case 10 %B-Spline function from patient data
%         numBreaths = simTime/((length(pmus_model)-1)/100);
%         lengthBreath = length(PmusHatSmoothNorm)-1;
%         startPos = 1;
%         for idx=1:simTime/4
%             Pmus(startPos:startPos + lengthBreath) = PmusHatSmoothNorm;
%             startPos = startPos + lengthBreath;
%             Pmus(startPos+1:startPos+231) = 0;
%             startPos = startPos + 231;
%         end
% end


