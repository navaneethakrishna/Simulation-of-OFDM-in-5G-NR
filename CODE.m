%% Description

clear;
close all;
clc;
%% Get the bandwidth, SCS, QAM and number of  RBs configuration information
[scs,bw] = input();                                             % get values of SCS and BW
Nrb_mat = [25,52,79,106,133,160,216,270,inf,inf,inf,...         % the resource block table for SCS and BW pairs
    inf;11,24,38,51,65,78,106,133,162,217,245,273;...
    inf,11,18,24,31,38,51,65,79,107,121,135];
bw_range = 1e6*[5,10,15,20,25,30,40,50,60,80,90,100];           % Bandwidths suppored by NR
scs_range = 1e3*[15,30,60];                                     % SCS options in this simulation
if (ismember(scs,scs_range) && ismember(bw,bw_range))
    Nrb = Nrb_mat((scs_range == scs),(bw_range == bw));         % get the number of resource blocks
else
    disp('Please enter valid values of SCS and BW');
    choice = input('Would you like to start over?(y/n)');
    if ((choice == 'y') || (choice == 'Y'))
        [scs,bw] = input();
    else
        return;    
    end      
end

if(Nrb == inf)                                                  % SCS and BW pair not defined
    str2 = 'Please choose a combination of SCS and BW for which the number of resource blocks are defined';
    fprintf(str2);
    choice = input('Would you like to start over?(y/n)');
    if ((choice == 'y') || (choice == 'Y'))
        [scs,bw] = input();
    else
        return;    
    end
end

%% Predefine parameters
Nfft = 4096;                                                    % Number of FFT taps
Ts = 1/(15e3*2048);                                             % LTE sampling time
Tc = 1/(scs*Nfft);                                              % NR sampling time
kappa = Ts/Tc;
mu = log2(kappa) - 1;
close figure 1;
prob = 1e-4;                                                    % Probabilty at whhich the value of PAPR is to be calculated from the CDF
band_used = Nrb*12*scs;                                         % Bandwidth used
no_frames = 100;                                                % Number of frames of data being transmitted
no_sub_car = Nrb*12;                                            % Number of subcarriers available for the chosen SCS and BW pair
no_slots = power(2,mu)*10;                                      % Number of slots per sub-frame
pow_tti = zeros(3,14*no_slots*10*no_frames);                    % total power transmitted

%% Obtain QAM scheme
prompt1 = ['Please enter the QAM scheme desired', newline, '2 : QPSK',newline,'16 : 16-QAM', newline,'64 : 64-QAM',newline,'all: to see a comparison between the three'];
dlgtit = 'QAM Scheme selection';
dims = [1 70];
config2 = inputdlg(prompt1,dlgtit,dims);
if (string(config2) == "all")
    modulation = [2,16,64];
elseif (str2double(config2) == 2) || (str2double(config2) == 16) || (str2double(config2) == 64) 
    modulation = str2double(config2);
else
    disp('Please enter a valid input in the next try! Bye bye!');
    return;
end

% Preallocate for speed
papr = zeros(length(modulation),14*no_slots*10*no_frames);      % matrix containing PAPR values of all OFDM samples
%% Full TX pipeline
for mod = 1 : length(modulation)
i = 1;
for frame = 1 : no_frames
    for sub_frame = 1 : 10
        for slot = 1 : no_slots
            for ofdm_sym_no = 1 : 14
                no_bits = no_sub_car*log2(modulation(mod));
%% Generating modulated symbols and OFDM signals

                seq = randi([0,1],1,no_bits);
                mod_sym = zeros(1,no_sub_car);
                if (modulation(mod) == 2)
                    for j = 0:(no_bits/modulation(mod)) -1
                        mod_sym(j+1) = ((1 - 2*seq(2*j+1)) + (1 - 2*seq(2*(j+1)))*1i)/sqrt(2);
                    end
                elseif (modulation(mod) == 16)
                    for j = 0:(no_bits/modulation(mod)) -1
                        mod_sym(j+1) = ((1-2*seq(4*j+1))*(2-(1-2*seq(4*j+1 + 2)))+((1-2*seq(4*j+1+1))*(2 - (1-2*seq(4*j+1+3))))*1i)/sqrt(10);
                    end
                else
                    for j = 0:(no_bits/modulation(mod)) -1
                        mod_sym(j+1) = (((1 - 2*seq(6*j + 1))*(4 - (1 - 2*seq(6*j + 3))*(2 - (1 - 2*seq(6*j + 5))))) + ((1 - 2*seq(6*j + 2))*(4 - (1 - 2*seq(6*j+4))*(2 - (1 - 2*seq(6*j + 6)))))*1i)/sqrt(42);
                    end
                end
%% Placing the data samples at the middle of the IFFT input sequence and obtaining OFDM samples

                ifft_input = zeros(Nfft,1);
                if(rem(no_sub_car,2) == 0)
                    data_start = 2049 - (no_sub_car/2);
                    ifft_input(data_start : data_start + length(mod_sym) -1) = mod_sym;
                else
                    data_start = 2048 - floor(no_sub_car/2);
                    ifft_input(data_start : data_start + length(mod_sym) -1) = mod_sym;
                end
                ofdm_samples = ifft(ifft_input)';
                
%% Cyclic prefix calculation
                if((ofdm_sym_no == 1) || (ofdm_sym_no == 8))
                    cp_length = 288 + 16*kappa; %%reset, mistake in REPORT.pdf
                else
                    cp_length = 288; %%reset, mistake in REPORT.pdf
                end
                cp = ofdm_samples(1,end-cp_length + 1:end);
                pp = max(power(abs([cp,ofdm_samples]),2));      % peak power
                ap = sum(power(abs([cp,ofdm_samples]),2))/4096; % average power
                   
                to_be_added = [cp,ofdm_samples];                % OFDM symbol with CP - for future use in any use case
                pow_tti(mod,i) = sum(abs(to_be_added));    % avearage power in the symbol with CP
                papr(mod,i) = pp/ap;
                i = i + 1;              
            end
        end
    end
end
end

%% Finding the CCDF
papr_dB = 10*log10(papr);
mean_dB = mean(papr_dB,2);
var_dB = mean(power(papr_dB,2),2) - power(mean_dB,2);
stdD_dB = sqrt(var_dB);
finish = max(papr_dB,[],2);
start = min(papr_dB,[],2) - 0.1;
if (length(modulation) == 1)
    x1 = start(1,1) : 0.05 : finish(1,1);
    ccdf_dB1 = zeros(1,length(x1));
    for i = 1 : length(x1)
        ccdf_dB1(1,i) = numel(find(papr_dB(1,:) > x1(i)));
    end
    ccdf_dB1 = ccdf_dB1 / max(ccdf_dB1);
    [~,I1] = min(abs(ccdf_dB1- prob));
    figure;
    semilogy(x1,ccdf_dB1);
    xlabel('PAPR');
    ylabel('Probabilities');
    title(['CCDF of PAPR for ',num2str(modulation),'-QAM Modulation']);
    at_prob = x1(I1);
    tot_avg_pow = sum(pow_tti,2)/(no_slots*10*no_frames);
    % Printing output
    disp([newline 'Sampling time = ',num2str(Tc),'s',' Bandwidth used = ', num2str(band_used/1e6),'MHz',' Nrb = ',num2str(Nrb)]);
    disp(['Sampling rate = ',num2str(1/(Tc*1e6)),'MHz',' kappa = ',num2str(kappa)])
    disp(['Average power transmitted per TTI = ',num2str(tot_avg_pow(1,1)),' CCDF at 0.001 probability = ', num2str(at_prob),'Mean = ', num2str(mean_dB),'dB, Variance = ',num2str(var_dB),'dB']);
else
x1 = start(1,1) : 0.05 : finish(1,1);
x2 = start(2,1) : 0.05 : finish(2,1);
x3 = start(3,1) : 0.05 : finish(3,1);
ccdf_dB1 = zeros(1,length(x1));
ccdf_dB2 = zeros(1,length(x2));
ccdf_dB3 = zeros(1,length(x3));
for i = 1 : length(x1)
    ccdf_dB1(1,i) = numel(find(papr_dB(1,:) > x1(i)));
end
for i = 1 : length(x2)
    ccdf_dB2(1,i) = numel(find(papr_dB(2,:) > x2(i)));
end
for i = 1 : length(x3)
    ccdf_dB3(1,i) = numel(find(papr_dB(3,:) > x3(i)));
end
ccdf_dB1 = ccdf_dB1 / max(ccdf_dB1);
ccdf_dB2 = ccdf_dB2 / max(ccdf_dB2);
ccdf_dB3 = ccdf_dB3 / max(ccdf_dB3);
[~,I1] = min(abs(ccdf_dB1- prob));
[~,I2] = min(abs(ccdf_dB2- prob));
[~,I3] = min(abs(ccdf_dB3- prob));
at_prob = [x1(I1), x2(I2), x3(I3)];
figure;
semilogy(x1,ccdf_dB1);
hold on;
semilogy(x2,ccdf_dB2);
hold on;
semilogy(x3,ccdf_dB3);
legend('QPSK','16 - QAM','64 - QAM');
xlabel('PAPR(dB)');
ylabel('Probability');
title('CCDF of OFDM System for different modulation schemes');
% Printing output
tot_avg_pow = sum(pow_tti,2)/no_slots*10*no_frames;
disp([newline 'Sampling time = ',num2str(Tc),'s',' Bandwidth used = ', num2str(band_used/1e6),'MHz',' Nrb = ',num2str(Nrb)]);
disp(['Sampling rate = ',num2str(1/(Tc*1e6)),'MHz',' kappa = ',num2str(kappa)])
disp(['For QPSK: Avg Tx power per TTI= ',num2str(tot_avg_pow(1,1)),' PAPR at 0.0001 prob = ', num2str(at_prob(1)),'dB,',' Mean = ', num2str(mean_dB(1)),'dB, Variance = ',num2str(var_dB(1)),'dB']);
disp(['For 16 - QAM: Avg Tx power per TTI = ',num2str(tot_avg_pow(2,1)),' PAPR at 0.0001 prob = ', num2str(at_prob(2)),'dB,',' Mean = ', num2str(mean_dB(2)),'dB, Variance = ',num2str(var_dB(2)),'dB']);
disp(['For 64 - QAM: Avg Tx power per TTI= ',num2str(tot_avg_pow(3,1)),' PAPR at 0.0001 prob = ', num2str(at_prob(3)),'dB,',' Mean = ', num2str(mean_dB(3)),'dB, Variance = ',num2str(var_dB(3)),'dB']);
end
%% Function definitions

function [scs,bw] = input() % creates dialogue box to get the SCS and BW
str1 = ' The table of number of resource blocks for different combinations of subcarrier spacing (rows) and \n signal bandwidth (columns) is will be shown on hitting any key. Please enter \n the desired values in the dialogue box.\n';
fprintf(str1);
pause;
figure('Position',[ 300 , 2000 , 1400 , 640]);
imshow('Table1.png');
prompt = {'Enter the SCS (KHz)','Enter the BW(MHz)'};
dlgtitle = 'Configuration definition';
dims = [1 70];
config = inputdlg(prompt,dlgtitle,dims);
scs = 1e3*str2double(config(1));
bw = 1e6*str2double(config(2));
end  
