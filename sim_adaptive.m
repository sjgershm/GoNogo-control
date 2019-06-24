function data = sim_adaptive(param,stim,R)
    
    % Simulationg of adaptive Pavlovian-instrumental Go/NoGo model.
    %
    % USAGE: data = sim_adaptive(param,stim,R)
    %
    % INPUTS:
    %   param - vector encoding the model parameters:
    %           param(1) = inverse temperature
    %           param(2) = prior mean, instrumental controller
    %           param(3) = prior confidence, instrumental controller
    %           param(4) = prior mean, Pavlovian controller
    %           param(5) = prior confidence, Pavlovian controller
    %           param(6) = initial Pavlovian weight
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
    mq = param(2);   % prior mean, instrumental
    pq = param(3);   % prior confidence, instrumental
    mv = param(4);   % prior mean, Pavlovian
    pv = param(5);   % prior confidence, Pavlovian
    
    if length(param)>5
        w0 = param(6);  % initial Pavlovian weight
    else
        w0 = 0.5;
    end
    
    u = unique(stim);
    S = length(u);
    v = zeros(S,1) + mv;
    q = zeros(S,2) + mq;
    Mv = zeros(S,1) + pv;
    Mq = zeros(S,2) + pq;
    L = log(w0) - log(1-w0);
    [~,opt] = max(R,[],2);
    
    for n = 1:length(stim)
                
        s = stim(n);  % stimulus
        w = 1./(1+exp(-L));
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
        
        % update model posterior
        if r == 1
            L = L + log(v(s)) - log(q(s,a));
        else
            L = L + log(1-v(s)) - log(1-q(s,a));
        end
        
        % update reward predictions
        Mv(s) = Mv(s) + 1;
        Mq(s,a) = Mq(s,a) + 1;
        v(s) = v(s) + (r-v(s))/Mv(s);
        q(s,a) = q(s,a) + (r-q(s,a))/Mq(s,a);
        
    end