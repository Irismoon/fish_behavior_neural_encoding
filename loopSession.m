function loopSession(func)
if ~isa(func,'function_handle')
    error('input must be a function handle');
    return;
end
[sessionID,fishID] = getfish();
for iS = 1:length(sessionID)
    func(sessionID{iS},fishID{iS});
end
end

