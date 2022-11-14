function tx_packet =packet_transmission(wifi_data_packet,~)
[~,len] = size(wifi_data_packet);
%% magnitude error, phase error, freq offset and awgn noise. 
idle_samples= zeros(1,100);
wifi_data_packet = [idle_samples wifi_data_packet];
wifi_data_packet_with_mag_err = (10^-5)* wifi_data_packet;
wifi_data_packet_with_phase_err = exp(-1i*3*pi/4).*wifi_data_packet_with_mag_err;

%wifi_data_packet_with_phase_err = wifi_data_packet_with_mag_err;

for ctr=1:len+100
    wifi_data_packet_with_freq_err(ctr)= exp(-1j*2*pi*0.00017*ctr)*wifi_data_packet_with_phase_err(ctr);
end

%% generating random numbers with mean 0 and std deviation of 10^-22.

mn = 0;
std_deviation = 10^-11;
wg_noise = std_deviation.*randn(len+100,1) + mn;

wifi_data_packet_with_awgn = wifi_data_packet_with_freq_err+ transpose(wg_noise);
tx_packet= wifi_data_packet_with_awgn;
%% plotting spectrum of sts. The length of STS is 224 samples. with 100 idle samples,
% the spectrum needs to be plotted is between 101 and 101+224= 324
plot(power(abs(tx_packet(101:260)),2));
end
