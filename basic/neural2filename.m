function filename = neural2filename(queryfile)
if ismember(queryfile,{'CalTrace3_original'})
    filename = 'CalTrace3_original';
elseif ismember(queryfile,{'Coherence3','A3','L_temp3'})
    filename = 'Coherence3';
elseif ismember(queryfile,{'denoised_Ca_OASIS'})
    filename = 'denoised_Ca_OASIS';
elseif ismember(queryfile,{'DFoF'})
    filename = 'DFoF';
elseif ismember(queryfile,{'Firing_Rate'})
    filename = 'Firing_Rate';
elseif ismember(queryfile,{'spikes_OASIS'})
    filename = 'spikes_OASIS';
else 
    error('no proper query file found');
end
end

