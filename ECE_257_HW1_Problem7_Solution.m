clear all
%Global Variables
number_bits = 4176; % 4160;
pilots_sy = [7 ,21 ,44 ,58];
pilots = [complex(1,0) complex(1,0) complex(1,0) complex(1,0)];
nulls = [1 28:38];
data_sy = [2:6 8:20 22:27 39:43 45:57,59:64];

N_sym = ceil(number_bits/length(data_sy));

%Generate Data
random_bits = randi([0 1],1,number_bits);
disp(["Randomness of Bits: ", mean(random_bits==0)])
save_signal(random_bits, "/0/bit_stream");
% 0 -- Make random bit stream the same -------------------------------------------------

data = zeros(1,N_sym*length(data_sy));
data(find(random_bits == 1)) = complex(1,0);
data(find(random_bits == 0)) = complex(-1,0);
save_signal(data, "/1/a");
% 1.a -- BPSK Narrowband Modulation ----------------------------------------------------

data_mat = zeros(N_sym,64);
data_mat(:,data_sy) = reshape(data,N_sym,[]);
data_mat(:,pilots_sy) = repmat(pilots,N_sym,1);
for i = 1: size (data_mat ,1)
    data_mat(i ,:) = ifft(data_mat(i ,:) ,64);
end
data_mat = [data_mat(:,end-15:end) data_mat];
save_signal(data_mat, "/1/b");
% 1.b -- OFDM Wideband Modulation ------------------------------------------------------

%STF Sequence
sts_f = zeros(1,64);
sts_f(1:27) = [0 0 0 0 -1-1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0 1+1i 0 0 0 1+1i 0 0 0 1+1i 0 0];
sts_f(39:64) = [0 0 1+1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0 -1-1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0];
sts_t = ifft(sqrt(13/6).*sts_f , 64);
sts_t = sts_t(1:16);

%LTF sequence
lts_f = [0 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1]; 
lts_t = ifft(lts_f , 64);

%Packet
packet = [repmat(sts_t, 1, 10) lts_t(33:64) lts_t lts_t];
save_signal(packet, "/1/c");
% 1.c -- Packetized data with STF/LTF --------------------------------------------------
figure(1)
plot(abs(sts_t))
hold on 
title ('STF: Time Domain ')

packet = [zeros(1,100) packet reshape(transpose(data_mat),1,[])];
[pxx, w] = pwelch(packet, rectwin(64), [], 64);
figure(2)
plot(w/pi, 10*log10(pxx), '-o')
hold on
title("Full Signal: Frequency Domain")
% 1.d -- Plot STF Time Domain + Signal Frequency Domain (w/ 64 frequency bins) ---------
%% 

%Distorted Signal
step = 1:length(packet);
packet = packet.* 10^(-5);
save_signal(packet, "/2/1");
% 2.1 -- Magnitude Distortion ----------------------------------------------------------

packet = packet.* exp(-1i*3*pi/4);
save_signal(packet, "/2/2");
% 2.2 -- Phase Offset ------------------------------------------------------------------

packet = packet.* exp(-1i*2*pi*0.00017.*step);
save_signal(packet, "/2/3");
% 2.3 -- Frequency Offset --------------------------------------------------------------

% Provided variance of 1e-11 which is sigma squared - but that is too much noise
% Not taking square root means practically no noise (50dB SNR) but is the
% most logical alternative from the assignment instructions
sigma = 10^(-11); 
noise = sigma/2 * randn(1, length(packet)) + 1i * sigma/2 * randn(1, length(packet));
snr = 10*log10(abs(packet) / abs(noise));
disp(snr);
packet = packet + noise;
save_signal(packet, "/2/4");
% 2.4 -- AWGN --------------------------------------------------------------------------
figure(3)
plot(abs(packet))
hold on
title ('After Channel: Time Domain')
xlim([0 500])
% 2.5 -- Plot STF Time Domain After Channel Distortion ---------------------------------

%Packet Detection
Rx = zeros(1,length(packet)); 
Ex = zeros(1,length(packet));
for i=1:length(packet)-31
    Rx(i) = sum(packet(i+[0:15]).*conj(packet(i+[16:31]))); 
    Ex(i) = sum(packet(i+[0:15]).*conj(packet(i+[0:15])));
end
figure(4)
plot(abs(Ex))
hold on
plot(abs(Rx)) 
legend('Energy','Self-correlation')
title ('Self-Correlation ')
xlim([0 500])
save_signal(Ex, "/3/energy");
% 3 -- Energy --------------------------------------------------------------------------
save_signal(Rx, "/3/self_correlation");
% 3 -- Self Correlation ----------------------------------------------------------------
% 3 -- Plot Time Domain Self Correlation + Energy --------------------------------------
% 3 -- Write Index of Packet Start -----------------------------------------------------

cross = xcorr(packet,sts_t);
cross = cross(length(packet):end); 
ind=find (abs(cross)>0.9*max(abs(cross))) ; 
figure(5)
plot(abs(cross))
title ('Cross-Correlation ');
xlim([0 500]);
save_signal(cross, "/4/cross_correlation");
% 4 -- Cross Correlation ---------------------------------------------------------------
% 4 -- Plot Cross Correlation in Time Domain -------------------------------------------

%Frequency Offset
lts_start = ind(1) + 160;
lts_1 = packet(lts_start:lts_start+63);
lts_2 = packet(lts_start+64:lts_start+127);
f_off = sum(imag(lts_2./lts_1))/(2*pi*64*64);
packet = packet.*exp(-1i*2*pi*f_off .*step); 
disp(['frequency offset estimation:',num2str(-f_off)])
save_signal(packet, "/5/a");
% 5a -- Frequency Offset ---------------------------------------------------------------

%Channel Estimation
lts_1 = packet(lts_start+32:lts_start+95);
lts_2 = packet(lts_start+96:lts_start+159);
lts_1_fft = fft(lts_1 ,64);
lts_2_fft = fft(lts_2 ,64);
H = (lts_1_fft+lts_2_fft)/2.*conj(lts_f);
disp (['magnitude estimation : ',num2str(abs(mean(H/(lts_f.* lts_f ))))])
save_signal(H, "/5/b");
% 5b -- Channel Estimation (Phase/Magnitude Distortion to Each Subcarrier) -------------

%Decoding
packet = packet(ind(1)+320:end); 
packet_mat = reshape(packet ,80 ,[]) ;
packet_mat = transpose(packet_mat);
packet_mat = packet_mat(:,17:end);

for i =1:size(packet_mat,1)
    packet_mat(i ,:) = fft(packet_mat(i ,:) ,64) ; 
end

packet_mat = packet_mat ./H;
packet_data = packet_mat(:,data_sy);
packet_data = reshape(packet_data ,1 ,[]) ; 
packet_ans = sign(real(packet_data));
rx_bits = packet_ans > 0;
save_signal(rx_bits, "/5/c")
% 5c -- Recovered Bits -----------------------------------------------------------------

count = length(find(rx_bits == random_bits));

ber = 1-count/(length(packet_data));
disp(['Bit error rate:',num2str(ber)]) 
figure(6)
scatter(real(packet_data) ,imag(packet_data)) 
title ('BPSK IQ')
xlim([-1.2, 1.2])
ylim([-1.2 1.2])

% -------------------------------------------------------------------------------------
function save_signal(s, name)
   % Save as real/complex in an array
   % IQ, T
   % Flattening so its always time
   s = reshape(s.', 1, []);
   out = [real(s); imag(s)];
   disp([name, size(s), size(out)])
   h5create("gold.h5", name, size(out));
   h5write("gold.h5", name, out);
end

function s = read_signal(name)
    out = h5read("gold.h5", name);
    s = out(1, :) + 1i * out(2, :);
end
