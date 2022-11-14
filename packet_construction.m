function [bit_stream,bpsk_mod_bit_stream,ofdm_sym_mat,ofdm_mod_sig,wifi_packet,preamble_seq]= packet_construction(LENGTH)
%% geneerating random bits of length 4160 bits
bit_stream=randn(1,LENGTH+16);
for ctr=1:LENGTH+16
    if(bit_stream(ctr)<0)
        bit_stream(ctr)=0;
    else
        bit_stream(ctr)=1;
    end
end

%% BPSK modulation 

bpsk_mod_bit_stream = (1== bit_stream).*1 + (0== bit_stream).*-1;

%% scatterplot(bpsk_mod_bit_stream);

%% dividing the symbols stream and then perform 64 point IFFT 
%  Dividing 4160+16 symbols by 52( N_data_subcarriers + N_pilot_subcarriers), we get 80 OFDM symbols.
% Nsd=48; Nsym=80. 
% d(k,n)= d(k+Nsd*n) 
    % k= 1....(Nsd); n= 1,2...,(Nsym)
Nsd= 48;
N_sym = (LENGTH+16)/Nsd; 
LENGTH=LENGTH+16;

% converting seriel data into parallel data
ofdm_sym = reshape(bpsk_mod_bit_stream,[Nsd,N_sym]);

%% Puncturing all data points at indices 7,21,35,49 and replacing it to be as +1+0j pilots
ofdm_sym_mat = zeros(64,N_sym);
ofdm_sym_mat(8,:) = 1+0i;
ofdm_sym_mat(22,:) = 1+0i;
ofdm_sym_mat(44,:) = 1+0i;
ofdm_sym_mat(58,:) = 1+0i;

ofdm_sym_mat(2:7,:) = ofdm_sym(1:6,:);
ofdm_sym_mat(9:21,:) = ofdm_sym(7:19,:);
ofdm_sym_mat(23:27,:) = ofdm_sym(20:24,:);
ofdm_sym_mat(39:43,:) = ofdm_sym(25:29,:);
ofdm_sym_mat(45:57,:)=ofdm_sym(30:42,:);
ofdm_sym_mat(59:64,:)=ofdm_sym(43:48,:);

%% Perform 64 point IFFT to OFDM modulate the signal. The signal is in Time domain now. 

ofdm_mod_sig = zeros(64,N_sym);
for i=1:N_sym
    ofdm_mod_sig(:,i) = ifft(ofdm_sym_mat(:,i),64);
end 

%FFT magnitude plot. 

for i=1:N_sym
    ofdm_spectrum(:,i)= fft(ofdm_mod_sig(:,i),64);
    plot(abs(ofdm_spectrum));
end 


%% Adding cyclic prefix to the ofdm modulated symbols. 
% Length of ofdm symbol is 80 : Data ( 64 ) + CP (16)
ofdm_with_cp = zeros(80,87);
ofdm_with_cp(1:16,:) = ofdm_mod_sig(49:64,:);
ofdm_with_cp(17:80,:) = ofdm_mod_sig;
%% pramble addition to the packet constructed. 
preamble_seq= preamble();
ofdm_data = reshape(ofdm_with_cp,1,[]);
wifi_packet = horzcat(preamble_seq,ofdm_data);

end









