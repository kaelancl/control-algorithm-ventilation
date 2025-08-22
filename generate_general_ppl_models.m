
%clear all
%close all
addpath("C:\Users\292636H\OneDrive - Curtin\PhD\Matlab")
addpath("C:\Users\292636H\OneDrive - Curtin\PhD\Matlab\B-MIQP Working Folder\.mat files")
%% Extract Single Breath\load patient data
load('patient_data/20_breaths_pat_5_1607_0004.mat')

R = 16.25;
E = 36.56;
PEEP = 8.5;
start_insp = [209
    645
    1226
    1767
    2243
    2702
    3247
    3745
    4282
    4756
    5241
    5715
    6206
    6728
    7203
    7753
    8258
    8786
    9257
    9737
    10222
    10699
];

%% Determine start of each patient inspiration
k_soe = find(diff(dataSignals.phase) == 32)+1;
Ecw = 19.56; %chest wall elastance
num_breaths = 20;
dataSignals.PEEP = 8.5;
%Determine start of each breath
% [pesShiftedEachBreath,pes_filt_shift, pes_shift,Pmus_inv_allBreaths, ...
%     volume_allBreaths, flow_allBreaths, paw_allBreaths, Pmus_inv_raw_allBreaths, phase_allBreaths] ...
%     = startOfInspPatFun(dataExtract,Ecw,num_breaths);

% blank arrays
coeff_bs_norm_array = zeros(20,30);
coeff_bs_array = zeros(20,30);
length_vent_insp_array = zeros(20,1);
insp_length_array = zeros(20,1);
Ti_array = zeros(20,1);
Paw_peak_array = zeros(20,1);
peaks_array = zeros(20,1);
peak_time_array = zeros(20,1);
counter =  1;
i = 1;
X_array = [];
params_array = [];
%beta_array = [];
%collect 20 breaths per patient
tiledlayout(5,4)
nexttile
while counter < num_breaths +1

    %if i is less than total number of breaths
    %if i <= length(volume_allBreaths)
        % volume = volume_allBreaths{i};
        % flow = flow_allBreaths{i};
        % paw = paw_allBreaths{i};
        % phase = phase_allBreaths{i};
   for i = 1:20
        % volume = dataSignals.vol(start_insp(i):start_insp(i+1));
        % flow = dataSignals.flow(start_insp(i):start_insp(i+1));
        % paw = dataSignals.pres(start_insp(i):start_insp(i+1));
        % phase = dataSignals.phase(start_insp(i):start_insp(i+1));
        
        volume = dataSignals.vol(start_insp(i):k_soe(i))/1000;
        flow = dataSignals.flow(start_insp(i):k_soe(i))/10000;
        paw = dataSignals.pres(start_insp(i):k_soe(i));
        phase = dataSignals.phase(start_insp(i):k_soe(i));

        %if no end expiratory volume present
        if volume(1) < 500
            %pes = pesShiftedEachBreath{i};
            %pmus = Pmus_inv_allBreaths{i};
            %pes = dataSignals.pes(start_insp(i):start_insp(i+1));
            %pes = pes - pes(1);
            pmus = paw - R*flow-E*volume-PEEP;
            % Peaks + lengths
           
            if length(find(diff(phase) == 32)) == 1
                [peaks_array(counter,1),peak_time_array(counter,1)] = min(pmus);
                length_vent_insp_array(counter,1) = find(diff(phase) == 32)+1;
                insp_length_array(counter,1) = find(pes<0,1,'last');
    
                if insp_length_array(counter,1) > length_vent_insp_array(counter,1)
                    pmus_vent = pmus(1:length_vent_insp_array(counter,1));
                    insp_length_array(counter,1) = find(pes_vent < 0,1,'last');
                end
    
                %normalise muscle pressure
                peak_pmus = min(pmus);
                pmus = pmus/peak_pmus;
                %peaks_array(i) = peak_pmus;
                %length of breath
                length_breath = length(flow)/100;
                %length_breath = length_vent_insp_array(counter,1)/100;
                totalTime = [0.01:0.01:length_breath]; %create time vector equal to length of inspiratory effort
    
                %create B-Spline model
                Kw = 0.25; %knot width, s
                Tmax = Kw * ceil(length_breath/Kw); %Time span of B-spline functions
                d = 3; %degree of B-splines
                M = Tmax/Kw + d; %number of B-splines
                %M = 10;
                knots = 0:Kw:Tmax; %knot vector
                knots = [0 0 0 knots Tmax Tmax Tmax]; %knot vector adding dummy values at start and end
    
                %create 2nd order B-spline fucnctions
                phi_2 = zeros(cast(M,"int32"),length(totalTime));
                for idx1 = 1:M
                    phi_2(idx1,(1:length(totalTime))) = bspline_basis(idx1,3,knots,totalTime);
                end
    
                % Model Pmus using B-Splines
                %lb = ones(M,1)*-100;
                lb = zeros(M,1);
                coeff_bs = lsqlin(phi_2',pmus,[],[],[],[],lb,[]);
    
                pmus_est = coeff_bs'*phi_2;

                figure()
                plot(pmus)
                hold on
                plot(pmus_est)
                hold off
                %% Invasive Resp Mech Estimates
                % C = [flow, volume];
                % d = [paw - pes - dataSignals.PEEP];
                % params = lsqlin(C,d,[],[],[],[],lb);
                % params_array = [params_array params];
                % 
                % Normalise B-Spline coefficients
    
                coeff_bs_norm = coeff_bs/max(coeff_bs);
    
    
                if length(coeff_bs_norm) >= 30
                    coeff_bs_norm_array(counter,:) = coeff_bs_norm(1:30);
                    coeff_bs_array(counter,:) = coeff_bs(1:30);
                else
                    coeff_bs_norm_array(counter,:) = [coeff_bs_norm',zeros(1,max(0,30-length(coeff_bs_norm)))];
                    coeff_bs_array(counter,:) =[coeff_bs',zeros(1,max(0,30-length(coeff_bs)))];
                end
                counter = counter + 1;
            end
        i = i+1;
    % else
    %     break
    %     counter
        end
    end
end


figure()
boxplot(coeff_bs_array)
hold on
plot(mean(coeff_bs_array),'.')

% figure()
% boxplot(coeff_bs_norm_array)
% hold on
% plot(mean(coeff_bs_norm_array),'.')

%length of breath
length_breath = 7;
%length_breath = length_vent_insp_array(counter,1)/100;
totalTime = [0.01:0.01:length_breath]; %create time vector equal to length of inspiratory effort

%create B-Spline model=
Kw = 0.25; %knot width, s
Tmax = Kw * ceil(length_breath/Kw); %Time span of B-spline functions
d = 3; %degree of B-splines
M = Tmax/Kw + d-1; %number of B-splines
%M = 10;
knots = 0:Kw:Tmax; %knot vector
knots = [0 0 0 knots Tmax Tmax Tmax]; %knot vector adding dummy values at start and end

phi_2 = zeros(cast(M,"int32"),length(totalTime));
for idx1 = 1:M
    phi_2(idx1,(1:length(totalTime))) = bspline_basis(idx1,3,knots,totalTime);
end
% median value of each normalised b-spline
coeff_med_norm = median(coeff_bs_array);
%remove negative spline coefficients
ind_coeff_pos = find(coeff_med_norm > 0);
coeff_med_norm(~ind_coeff_pos) = 0;
%%force monotonicity
%%*******************************************************************

coeff_med_norm(1:3) = 0;
pmus_model = coeff_med_norm*phi_2;


pmus_data = paw-R*flow/1000-E*volume/1000-PEEP;
pmus_data_peak = min(pmus_data);
pmus_data = pmus_data/pmus_data_peak;

%normalise model
pmus_model = pmus_model/max(pmus_model);


figure()
plot(pmus_data)
hold on
plot(pmus_model)

figure(); %scatter(linspace(0,25*10,10),coeff_med_norm(1:10))


hold on; plot(pmus_model)
%plot splines
for i = 1:(size(coeff_med_norm,2))
    plot(coeff_med_norm(i)*phi_2(i,:));
end

%% Filter Ppl model signal
% Fs = 100;            % Sampling frequency
% dt = 1/Fs;             % Sampling period
% n = length(pes_model);             % Length of signal
% 
% %perform fourier transform
% fhat = fft(pes_model);
% 
% %determine frequencies
% freq = Fs/n*(0:n-1);
% %create system variables
% sys = frd(fhat,freq,'FrequencyUnit','Hz');
% 
% %fft plot
% figure()
% plot(Fs/n*(0:n-1),abs(fhat),"-")
% hold on

