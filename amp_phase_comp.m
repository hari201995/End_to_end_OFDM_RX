function X_cap= amp_phase_comp(rx_packet,H_k)
[row,col]=size(rx_packet);
X_cap=zeros(row,col);
%% amplitude compensation 
% multiply every ofdm subcarrier with H_cap
for ctr=1:col
   temp=transpose(rx_packet(:,ctr))./H_k;
   X_cap(:,ctr) = temp;
end
inf_val_Arr=[1,28:38];
X_cap(inf_val_Arr,:)=0;
%scatterplot(X_cap(:,1:end));
end