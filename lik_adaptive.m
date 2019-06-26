function [lik, latents] = lik_adaptive(param,data)
    
    bt = param(1);  % inverse temperature
    
    if length(param) > 2
        mq = param(2);   % prior mean, instrumental
        pq = param(3);   % prior confidence, instrumental
        mv = param(4);   % prior mean, Pavlovian
        pv = param(5);   % prior confidence, Pavlovian
    else
        mq = 0.5; mv = 0.5;
        pq = 2; pv = 2;
    end
    
    if length(param)>6
        w0 = param(7);  % initial Pavlovian weight (log odds)
    elseif length(param)==3
        w0 = param(3);
    else
        w0 = 0.5;
    end
    
    u = unique(data.s);
    S = length(u);
    v = zeros(S,1) + mv;
    q = zeros(S,2) + mq;
    Mv = zeros(S,1) + pv;
    Mq = zeros(S,2) + pq;
    lik = 0;
    
    for n = 1:data.N
    
        if n == 1 %|| data.block(n)~=data.block(n-1)
            L = log(w0) - log(1-w0);
        end
        
        s = data.s(n);  % stimulus
        w = 1./(1+exp(-L));
        d = (1-w)*q(s,1) - (1-w)*q(s,2) - w*v(s);
        P = 1./(1+exp(-bt*d)); % probability of NoGo
        c = data.a(n);
        r = data.r(n);
        
        if c==1
            lik = lik + safelog(P);
        else
            lik = lik + safelog(1-P);
        end
        
        if nargout > 1
            latents.q(n,:) = q(s,:);
            latents.v(n,1) = v(s);
            latents.w(n,1) = w;
            latents.P(n,1) = P;
            if (data.acc(n)==1 && c==1) || (data.acc(n)==0 && c==2)
                latents.acc(n,1) = P;
            else
                latents.acc(n,1) = 1-P;
            end
        end
        
        if r == 1
            L = L + log(v(s)) - log(q(s,c));
        else
            L = L + log(1-v(s)) - log(1-q(s,c));
        end
        
        Mv(s) = Mv(s) + 1;
        Mq(s,c) = Mq(s,c) + 1;
        v(s) = v(s) + (r-v(s))/Mv(s);
        q(s,c) = q(s,c) + (r-q(s,c))/Mq(s,c);
        
    end