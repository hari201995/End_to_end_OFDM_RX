function [packet_index, packet_det]= packet_detection(tx_packet,~)
m=17;
%% self correlator
%% comparing successive energies 
% [~,len_tx_pkt] = size(tx_packet);
for ctr=2:length(tx_packet)-32
    R_sum = 0;
    E_sum = 0;
    for pivot=m:m+15
        R_sum = R_sum + tx_packet(pivot)*(tx_packet(pivot-16)');
        E_sum = E_sum + tx_packet(pivot)*(tx_packet(pivot)');
        %x_m(pivot)= tx_packet(1,pivot)*tx_packet(1,pivot-16)';
    end
    r_m(m-16) = R_sum;
    e_m(m-16) = E_sum;
     m=ctr+16;
end
plot(abs(r_m));
hold on
plot(abs(e_m));
hold off

%% from the plot, we see the power of R_m and E_m is same at 101st sample.
packet_index=101;
packet_det=1;

end