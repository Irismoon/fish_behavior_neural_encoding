function roi = pixel2region(x,y,z,sessionID,fishID)

prepath = fullfile(getpath('neural activity',sessionID,fishID),'Coherence3*');
fileinfo = dir(prepath);
[~,I] = max([fileinfo(:).datenum]);
load(fullfile(fileinfo(I).folder,fileinfo(I).name),'L_temp3');
roi = L_temp3(y,x,z);%note that since the x, y is reverse from pixel to region, i.e., in figure,
%[X,Y] is actually [column, row] of matrix, thus here y and x is reverse. 
end

