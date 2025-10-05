%% ===================== Wrapper: DEsyst with constant ca & G_syn over window =====================
function dy = NVC_chunk(time, state, ca_const, Gsyn_const, tmax, t_start, t_end)
% Reuse your NVC exactly, but pass *constant* signals as 1-element arrays.
% DEsyst internally indexes into these arrays; with length=1 it always selects the single value.
ca_arr   = ca_const;     % scalar; NVC treats this as length-1 array
G_arr    = Gsyn_const;   % scalar; NVC treats this as length-1 array
dy = NVC(time, state, ca_arr, G_arr, tmax, t_start, t_end);
end