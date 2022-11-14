function H_cap = channel_estimation(freq_comp_tx_packet,stf_end)

L_seq =[1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, 1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, 0, ... 
1, -1, -1, 1, 1, -1, 1, -1, 1, -1, -1, -1, -1, -1, 1, 1, -1, -1, 1, -1, 1, -1, 1, 1, 1, 1];
X_k = [L_seq(27),L_seq(28:53),zeros(1,11),L_seq(1:26)]; 

ltf_start=stf_end+1;
ltf_end = ltf_start+160-1;

L_seq_t_domain = freq_comp_tx_packet(ltf_start+32:ltf_end);
L_seq_sym1 = L_seq_t_domain(1:64);
%% FFT 
Y_k = fft(L_seq_sym1,64);
H_cap = Y_k./X_k;
inf_val_Arr=[1,28:38];
H_cap(inf_val_Arr)=0;
mag_spec=abs(H_cap);
plot(mag_spec);
phase_spec = angle(H_cap);
plot(phase_spec);
end