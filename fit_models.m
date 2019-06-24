function [results, bms_results] = fit_models(data,models,results)
    
    likfuns = {'lik_fixed' 'lik_adaptive' 'lik_adaptive_reduced'};
    
    if nargin < 1; data = load_data; end
    if nargin < 2; models = 1:length(likfuns); end
    
    pmin = 0.01;
    pmax = 80;
    btmax = 30;
    
    for mi = 1:length(models)
        m = models(mi);
        disp(['... fitting model ',num2str(m)]);
        fun = str2func(likfuns{m});
        
        switch likfuns{m}
            
            case 'lik_fixed'
                
                param(1) = struct('name','invtemp','logpdf',@(x) 0,'lb',1e-3,'ub',btmax);
                param(2) = struct('name','w','logpdf',@(x) 0,'lb',0.001,'ub',0.999);
                param(3) = struct('name','mq','logpdf',@(x) 0,'lb',0.001,'ub',0.999);
                param(4) = struct('name','pq','logpdf',@(x) 0,'lb',pmin,'ub',pmax);
                param(5) = struct('name','mq','logpdf',@(x) 0,'lb',0.001,'ub',0.999);
                param(6) = struct('name','pq','logpdf',@(x) 0,'lb',pmin,'ub',pmax);
                
            case 'lik_adaptive'
                
                param(1) = struct('name','invtemp','logpdf',@(x) 0,'lb',1e-3,'ub',btmax);
                param(2) = struct('name','mq','logpdf',@(x) 0,'lb',0.001,'ub',0.999);
                param(3) = struct('name','pq','logpdf',@(x) 0,'lb',pmin,'ub',pmax);
                param(4) = struct('name','mq','logpdf',@(x) 0,'lb',0.001,'ub',0.999);
                param(5) = struct('name','pq','logpdf',@(x) 0,'lb',pmin,'ub',pmax);
                param(6) = struct('name','w0','logpdf',@(x) 0,'lb',0.001,'ub',0.999);
                
            case 'lik_adaptive_reduced'
                
                param(1) = struct('name','invtemp','logpdf',@(x) 0,'lb',1e-3,'ub',btmax);
                param(2) = struct('name','mq','logpdf',@(x) 0,'lb',0.001,'ub',0.999);
                param(3) = struct('name','pq','logpdf',@(x) 0,'lb',pmin,'ub',pmax);
                param(4) = struct('name','mq','logpdf',@(x) 0,'lb',0.001,'ub',0.999);
                param(5) = struct('name','pq','logpdf',@(x) 0,'lb',pmin,'ub',pmax);
                fun = str2func('lik_adaptive');
                
        end
        
        results(m) = mfit_optimize(fun,param,data);
        clear param
    end
    
    % Bayesian model selection
    if nargout > 1
        bms_results = mfit_bms(results,1);
    end