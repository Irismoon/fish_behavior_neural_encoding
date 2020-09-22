function adjustImagContrast(filename)
%this file is deprecated since it's out of memory
iminfo = imfinfo(filename);
numFrame = length(iminfo);

lohi = zeros(numFrame,2);
parfor iFrame = 1:numFrame
    lohi(iFrame,:) = stretchlim(imread(filename,iFrame));
end
lohi = quantile(lohi,0.8);
[filepath,name,ext] = fileparts(filename);

fprintf('start writing...');
for iFrame = 1:numFrame
    if mod(iFrame,500)==0 disp(iFrame); end
    imageTmp = imadjust(imread(filename,iFrame),lohi);
    if iFrame==1
		imwrite(imageTmp,fullfile(filepath,[name '_contrastAdjusted' ext]));
	elseif iFrame>1
		imwrite(imageTmp,fullfile(filepath,[name '_contrastAdjusted' ext]),'WriteMode','append');
	end
end
disp('done!');
end