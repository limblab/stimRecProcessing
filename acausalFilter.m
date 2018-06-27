function [filteredData] = acausalFilter(data)
    % data is an MxN matrix, Filtering occurs along the first dimension
    
    
    filteredData = [];
    if(isempty(data))
        return;
    end
    numPad = 300;
    % make filter
    [b,a] = butter(6,[500]/(30000/2),'high');
    % pad data
    data = [repmat(mean(data(1:min(100,size(data,1)),:)),[numPad,1]);...
        data;...
        repmat(mean(data(end-min(99,size(data,1)-1):end,:)),[numPad,1])];
    % acausal filter
    filteredData = fliplr(filter(b,a,fliplr(data')')')';
    
%     [b,a] = butter(2,[7500]/(30000/2),'low');
%     f = fliplr(filter(b,a,fliplr(f')')')';

    % remove padding
    filteredData = filteredData(numPad+1:end-numPad,:);
end

