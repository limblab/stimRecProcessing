function [f] = acausalFilter(data)
    % data is an MxN matrix, M = data, N = number of artifacts (channels,
    % repeated stimuli, etc.). Data has to be along the first dimension
    f = data;
    if(isempty(data))
        return;
    end
    numPad = 300;
    % make filter
    [b,a] = butter(6,[500]/(30000/2),'high');
    % pad data
    data = [mean(data(1:min(100,size(data,1)),:))+zeros(numPad,size(data,2));...
        data;...
        mean(data(end-min(99,size(data,1)-1):end,:))+zeros(numPad,size(data,2))];
    % acausal filter
    f = fliplr(filter(b,a,fliplr(data')')')';
    
%     [b,a] = butter(2,[7500]/(30000/2),'low');
%     f = fliplr(filter(b,a,fliplr(f')')')';

    % remove padding
    f = f(numPad+1:end-numPad,:);
end

