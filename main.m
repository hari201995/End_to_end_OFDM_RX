clear all
close all
%% calling packet construction
LENGTH = 4160;
[bit_stream,bpsk_mod_bit_stream,ofdm_sym_mat,ofdm_mod_sig,wifi_packet,pre_seq] = packet_construction(LENGTH);
noise = input('run with noise - 1 or without noise - 0 \n');

if noise == 1

    %% calling packet transmission 
    tx_packet = packet_transmission(wifi_packet,LENGTH);
    %scatterplot(tx_packet);
    %% Receiver construction 
    %sliding correlator and packet detection
    [packet_no,packet_det]= packet_detection(tx_packet,pre_seq);
    
    % cross corr for time sync. 
    packet_idx = packet_sync(tx_packet,pre_seq);
    
    %freq offset est. 
     delta_f = freq_offset_est(tx_packet,packet_idx,packet_idx+160-1);
    
    %% frequency compensation at Rx. 
    
    idx = 1:1:length(tx_packet);
    freq_comp_vec = exp(1j*2*pi*idx*delta_f);
    scatterplot(tx_packet(1,1:end));
    freq_comp_tx_packet = exp(-1j*2*pi*idx*delta_f).*tx_packet;
    scatterplot(freq_comp_tx_packet);
    %% Channel estimation
    H_k = channel_estimation(freq_comp_tx_packet,packet_idx+160-1);

end

if noise ==0
    freq_comp_tx_packet = wifi_packet;
    packet_idx=1;
end

% Decoding data points at Rx
if noise ==1
    rx_processed_pkt = reshape(freq_comp_tx_packet(packet_idx+320:end),[80,87]);
end
if noise ==0
    rx_processed_pkt = reshape(wifi_packet(321:end),[80,87]);
end
%rx_processed_pkt = reshape(tx_packet,[80,87]);

rx_carrier_demod = zeros(64,87);
%% CP removal
temp_fft_ip= rx_processed_pkt(17:80,:);

%% FFT demod
for pivot=1:87
    rx_carrier_demod(:,pivot) = fft(temp_fft_ip(:,pivot),64); 
end

% amplitude and phase compensation
if noise == 1
    X_cap= amp_phase_comp(rx_carrier_demod,H_k);
end
%PSK demodulator
%rx_seriel_data = packet_deconstructor(X_cap);
temp=real(rx_carrier_demod);
temp2 = round(real(rx_carrier_demod));
%% packet deconstruct
rx_seriel_data = packet_deconstructor(rx_carrier_demod);

check = rx_seriel_data - bit_stream;
M = find(check);
success_rate = (4176-length(M))/4176;











