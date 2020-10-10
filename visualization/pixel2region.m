function roi = pixel2region(x,y,z)
load(fullfile(getpath('neural activity',sessionID,fishID),'Coherence3'),'L_temp3');
roi = L_temp3(x,y,z);
end

