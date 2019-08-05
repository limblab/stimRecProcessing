function [filteredData] = acausalFilter(data)
    % data is an 2-D matrix, filters along the first dimension. 
    
    filteredData = [];
    if(isempty(data))
        return;
    end
    
    if(size(data,1) == 1 && size(data,2) ~= 1)
        warning('incorrectly sized matrix, assuming transpose is desired');
        data = data';
    end
    
    numPad = 300;
    % make filter
%     [b,a] = butter(6,[500]/(30000/2),'high');
    [b,a] = butter(4,[500,2500]/(30000/2),'bandpass');
    % pad data
    data = [repmat(mean(data(1:min(15,size(data,1)),:)),[numPad,1]);...
        data;...
        repmat(mean(data(end-min(15,size(data,1)-1):end,:)),[numPad,1])];
    % acausal filter
    filteredData = flip(filter(b,a,flip(data,1)),1);
    
%     [b,a] = butter(2,[7500]/(30000/2),'low');
%     f = fliplr(filter(b,a,fliplr(f')')')';

    % remove padding
    filteredData = filteredData(numPad+1:end-numPad,:);
end

