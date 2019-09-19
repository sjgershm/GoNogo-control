function data = sim_fixed(param,stim,R)
    
    % Simulationg of adaptive Pavlovian-instrumental Go/NoGo model.
    %
    % USAGE: data = sim_adaptive(param,stim,R)
    %
    % INPUTS:
    %   param - vector encoding the model parameters:
    %           param(1) = inverse temperature
    %           param(2) = Pavlovian weight
    %           param(3) = prior mean, instrumental controller
    %           param(4) = prior confidence, instrumental controller
    %           param(5) = prior mean, Pavlovian controller
    %           param(6) = prior confidence, Pavlovian controller
    %   stim - [N x 1] vector containing the sequence of stimuli (indicated by integers)
    %   R - [S x A] matrix of stimulus-action reward probabilities
    %   
    % OUTPUTS:
    %   data - structure containing simulation results, with the following fields:
    %           .s - stimulus sequence
    %           .a - action sequence
    %           .r - reward sequence
    %           .acc - accuracy (1 = optimal action chosen)
    %           .w - Pavlovian weight sequence
    %
    % DEMOS:
    %   stim = ones(100,1);
    %   param = [5 0.5 2 0.5 2 0.5];
    %   R = [0.5 0.5]; % uncontrollable environment
    %   data = sim_adaptive(param,stim,R);
    %   plot(data.w);
    %   R = [0.9 0.1]; % controllable environment
    %   data = sim_adaptive(param,stim,R);
    %   hold on; plot(data.w,'-r');
    %   xlabel('Trial'); ylabel('Pavlovian weight');
    %
    % Sam Gershman, Jan 2019
    
    bt = param(1);   % inverse temperature
    w = param(2);   % Pavlovian weight
    
    if length(param) > 2
        mq = param(3);   % prior mean, instrumental
        pq = param(4);   % prior confidence, instrumental
        mv = param(5);   % prior mean, Pavlovian
        pv = param(6);   % prior confidence, Pavlovian
    else
        mq = 0.5; mv = 0.5;
        pq = 2; pv = 2;
    end
    
    u = unique(stim);
    S = length(u);
    v = zeros(S,1) + mv;
    q = zeros(S,2) + mq;
    Mv = zeros(S,1) + pv;
    Mq = zeros(S,2) + pq;
    [~,opt] = max(R,[],2);
    
    for n = 1:length(stim)
                
        s = stim(n);  % stimulus
        d = (1-w)*q(s,1) - (1-w)*q(s,2) - w*v(s);
        P = 1/(1+exp(-bt*d)); % probability of NoGo
        if rand < P
            a = 1;
        else
            a = 2;
        end
        
        % sample reward
        r = rand < R(s,a);
        acc = a==opt(s);
        
        % store data
        data.a(n,1) = a;
        data.r(n,1) = double(r);
        data.acc(n,1) = double(acc);
        data.w(n,1) = w;
        data.s(n,1) = s;
        
        % update reward predictions
        Mv(s) = Mv(s) + 1;
        Mq(s,a) = Mq(s,a) + 1;
        v(s) = v(s) + (r-v(s))/Mv(s);
        q(s,a) = q(s,a) + (r-q(s,a))/Mq(s,a);
        
    end