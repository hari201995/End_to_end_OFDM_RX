function rx_seriel_data = packet_deconstructor(X_cap)
[row,cols] = size(X_cap);
%% removing an additional phase. 
bpsk_demod_mat = zeros(row,cols);
for k=1:cols
    for m=1:row
        if (angle(X_cap(m,k))>=0)
            bpsk_demod_mat(m,k)= 1;
        else
            bpsk_demod_mat(m,k) =-1;
        end
    end
end

scatterplot(bpsk_demod_mat(:,1));
Rx_data = zeros(64,87);

for j=1:87
    for i=1:64
        if round(real(X_cap(i,j)))== -1
            Rx_data(i,j) = 0;
        end
        if round(real(X_cap(i,j))) == 1
            Rx_data(i,j) = 1;
        end
    end
end


rx_data_bits = zeros(48,87);
rx_data_bits(1:6,:) = Rx_data(2:7,:);
rx_data_bits(7:19,:) = Rx_data(9:21,:);
rx_data_bits(20:24,:) = Rx_data(23:27,:);
rx_data_bits(25:29,:) = Rx_data(39:43,:);
rx_data_bits(30:42,:) = Rx_data(45:57,:);
rx_data_bits(43:48,:) = Rx_data(59:64,:);
rx_seriel_data = reshape(rx_data_bits,[1,48*87]);

end