%% �ֶ�ȥƽ��ֵ
%input: window size, data
%output: ȥƽ�����ֵ
function output=remove_average(window, data)
    frameNum=length(data);
    output=data;
    for mm=1:fix(frameNum/window)
        output(window*(mm-1)+1:window*mm)=data(window*(mm-1)+1:window*mm)...
            -mean(data(window*(mm-1)+1:window*mm)); 
        if(mm==fix(frameNum/window))
            output(window*mm+1:end)=data(window*mm+1:end)-mean(data(window*mm+1:end));
        end
    end

end