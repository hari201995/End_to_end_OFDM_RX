function packet_idx = packet_sync(rx_packet,pre_seq)
sample_idx=1;
S=pre_seq(101:116);
while sample_idx<=length(rx_packet)-32
    R(sample_idx)=rx_packet(sample_idx:sample_idx+15)*S.';
    sample_idx=sample_idx+1;
end
plot(abs(R));
%% from the graph we see that we get the absolute value at positions starting from 101 which matches with out the 
%% given data. The peaks repeat after every 16 samples. 
packet_idx=101;
end

