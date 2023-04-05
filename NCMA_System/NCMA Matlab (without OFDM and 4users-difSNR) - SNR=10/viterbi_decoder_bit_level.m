function [uncoded_out] = viterbi_decoder_bit_level(RawBits,trellis,pkt_length)

% This viterbi decoder works in forward direction. That is, it computes the
% branch metrics for the next states. Other viterbi decoders look back to
% calculate the path metric to the current state.

% Note that in feed forward viterbi, we can go to a particular state with
% only 1 input. That is, if bit 0 information bit pushes the encoder to
% first state, then bit 1 will not go to the first state.
no_states=size(trellis.outputs,1);
out=10*ones(no_states,pkt_length);
svSoft_v=zeros(no_states,pkt_length);
pmSoft = Inf (no_states,1);
pmSoft(1) = 0;
metric=[];
for ii = 1:pkt_length
    pmSoft_n = Inf (no_states,1);
    svSoft = zeros (no_states,1);
    BMSoft = [RawBits(4*(ii-1)+1)+RawBits(4*(ii-1)+3), RawBits(4*(ii-1)+1)+RawBits(4*(ii-1)+4), RawBits(4*(ii-1)+2)+RawBits(4*(ii-1)+3), RawBits(4*(ii-1)+2)+RawBits(4*(ii-1)+4)];

  
    for i = 1:no_states % i is the present state

        pm_0 = pmSoft(i) + BMSoft(trellis.outputs(i,1)+1); 
        pm_1 = pmSoft(i) + BMSoft(trellis.outputs(i,2)+1);
        
        [pmSoft_n(trellis.nextStates(i,1)+1), indx_a]=min([pm_0,pmSoft_n(trellis.nextStates(i,1)+1)]);
        [pmSoft_n(trellis.nextStates(i,2)+1), indx_b]=min([pm_1,pmSoft_n(trellis.nextStates(i,2)+1)]);

        if indx_a==1
            svSoft(trellis.nextStates(i,1)+1,1)  = i; 
            out(trellis.nextStates(i,1)+1,ii) = 0;
        end
        if indx_b==1
            svSoft(trellis.nextStates(i,2)+1,1)  = i; 
            out(trellis.nextStates(i,2)+1,ii) = 1;       
        end
    end
    pmSoft=pmSoft_n;
    metric=[metric,pmSoft_n];
    svSoft_v(:,ii) = svSoft;
end

uncoded_out= 10*ones(pkt_length,1);
%[x, current_state] = min(pmSoft);
current_state=1; % Zero state is forced at tx side.
for ii= pkt_length:-1:1
    uncoded_out(ii,1) = out(current_state,ii);
    current_state = svSoft_v(current_state,ii);
end