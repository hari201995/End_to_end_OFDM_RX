function delta_f = freq_offset_est(tx_packet,stf_start,stf_end)
    ltf_start=stf_end+1;
    ltf_end = ltf_start+160-1;
    L_seq_t_domain = tx_packet(ltf_start+32:ltf_end);
    L_seq_sym1 = L_seq_t_domain(1:64);
    L_seq_sym2 = L_seq_t_domain(65:128);
    delta_f = L_seq_sym1./L_seq_sym2;
    delta_f = imag(delta_f);
    delta_f_deg = (1/(2*pi*64))*delta_f;
    delta_f = mean(delta_f_deg);
end