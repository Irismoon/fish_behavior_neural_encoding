% from the selected good bouts, decode the prey location in the visual field in a time-varying way.

pcoi = [1 2 3 4 5];
arrayfun(@(i) generate_MIP(W(:,i),csessionID,cfishID,['_' num2str(i) 'ndPC_random']),pcoi,'un',0);



%%
proj_focus

