function data = load_data(filename)
    
    % Load data from Go/NoGo experiments.
    %
    % USAGE: data = load_data(filename)
    %
    % INPUTS:
    %   filename - string specifying the csv file to be read
    %
    % OUTPUTS:
    %   data - [1 x nSubjects] array of stuctures, where each element
    %           corresponds to a single subject's data
    %
    % Sam Gershman, January 2019
    
    f = fopen(filename);
    y = fgetl(f);
    C = strsplit(y,',');
    fclose(f);
    X = csvread(filename,1);
    
    k = strcmp(C,'sub');
    S = unique(X(:,k));
    
    for s = 1:length(S)
        ix = X(:,k)==S(s);
        
        for i = 1:length(C)
            if i~=k
                data(s).(C{i}) = X(ix,i);
            end
        end
        
        data(s).N = length(data(s).a);
        acc(s) = nanmean(data(s).acc);
        
        for i = unique(data(s).s)'
            A(i) = mean(data(s).acc(data(s).s==i));
        end
        bad(s,1) = any(A<0.3);
    end
    
    data(acc'<0.5 | bad) = [];