function simdata = sim_bias_variance
    
    data = load_data('data1.csv');
    load results1
    
    P{1} = [0.25 0.75; 0.75 0.25; 0.5 0.5]; % low control
    P{2} = [0.25 0.75; 0.75 0.25; 0.2 0.8]; % high control
    
    for s = 1:length(data)
        c = data(s).cond(1);
        param = results(3).x(s,:);
        simdata(s,1) = sim_adaptive(param,data(1).s,P{c});
        param = results(1).x(s,:);
        param(2) = 0;
        simdata(s,2) = sim_fixed(param,data(1).s,P{c});   % pure instrumental
        param(2) = 1;
        simdata(s,3) = sim_fixed(param,data(1).s,P{c});   % pure Pavlovian
    end
    
    save simdata simdata