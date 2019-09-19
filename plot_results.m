function plot_results(fig,results)
    
    rng(1); % set random seed for reproducibile bootstrapped confidence intervals
    
    switch fig
        
        case 'weight_dynamics'
            
            load results1
            data = load_data('data1.csv');
            
            W = cell(1,2);
            for s = 1:length(data)
                c = data(s).cond(1);
                ix = data(s).s==3;
                W{c}(end+1,:) = results(3).latents(s).w(ix);
            end
            
            plot([W{1}(1,:)' W{2}(1,:)'],'LineWidth',4);
            set(gca,'FontSize',20,'YLim',[0 1]);
            ylabel('Adaptive weight (w)','FontSize',25);
            xlabel('Trial','FontSize',25);
            legend({'Low control' 'High control'},'FontSize',25,'Location','NorthOutside');
            
        case 'gobias_noexclude'
            
            acc = cell(2,2);
            for j = 1:2
                if j==1
                    data = load_data('data1.csv',0);
                else
                    data = load_data('data2.csv',0);
                end
                for s = 1:length(data)
                    c = data(s).cond(1);
                    go = data(s).s==1;
                    nogo = data(s).s==2;
                    acc{j,c}(end+1) = mean(data(s).acc(go)) - mean(data(s).acc(nogo));
                end
            end
            
            for c = 1:2
                for j = 1:size(acc,1)
                    m(c,j) = mean(acc{j,c});
                    err(c,j) = std(acc{j,c})./sqrt(length(acc{j,c}));
                end
            end
            
            barerrorbar(m',err');
            set(gca,'XTickLabel',{'Experiment 1' 'Experiment 2'},'FontSize',25,'XLim',[0.5 2.5],'YLim',[0 0.55]);
            ylabel('Go bias','FontSize',25);
            legend({'Low control' 'High control'},'FontSize',25);
            
        case 'gobias'
            
            subplot(2,2,1); plot_results('gobias1'); title('Experiment 1','FontSize',25,'FontWeight','Bold');
            subplot(2,2,2); plot_results('gobias2'); title('Experiment 2','FontSize',25,'FontWeight','Bold');
            subplot(2,2,3); plot_results('weight_gobias1');
            subplot(2,2,4); plot_results('weight_gobias2'); legend('off');
            set(gcf,'Position',[200 200 850 650]);
        
        case 'gobias1'
            
            % Go bias for experiment 1
            
            data = load_data('data1.csv');
            load results1
            
            acc = cell(2,2);
            for s = 1:length(data)
                c = data(s).cond(1);
                go = data(s).s==1;
                nogo = data(s).s==2;
                acc{1,c}(end+1) = mean(data(s).acc(go)) - mean(data(s).acc(nogo));
                acc{2,c}(end+1) = mean(results(3).latents(s).acc(go)) - mean(results(3).latents(s).acc(nogo));
            end
            
            for c = 1:2
                for j = 1:size(acc,1)
                    m(c,j) = mean(acc{j,c});
                    err(c,j) = std(acc{j,c})./sqrt(length(acc{j,c}));
                end
            end
            
            barerrorbar(m',err');
            set(gca,'XTickLabel',{'Data' 'Model'},'FontSize',25,'XLim',[0.5 2.5],'YLim',[0 0.55]);
            ylabel('Go bias','FontSize',25);
            legend({'Low control' 'High control'},'FontSize',25);
            
            [~,p,~,stat] = ttest2(acc{1,1},acc{1,2});
            d = (mean(acc{1,1}) - mean(acc{1,2}))./sqrt(0.5*(var(acc{1,1})+var(acc{1,2})));
            disp(['go bias (experiment 1): t(',num2str(stat.df),') = ',num2str(stat.tstat),', p = ',num2str(p),', d = ',num2str(d)]);
            
        case 'gobias2'
            
            % Go bias for experiment 2
            
            data = load_data('data2.csv');
            load results2
            
            for s = 1:length(data)
                acc(s,1) = mean(data(s).acc(data(s).s==1)) - mean(data(s).acc(data(s).s==2));
                acc(s,2) = mean(data(s).acc(data(s).s==4)) - mean(data(s).acc(data(s).s==5));
                macc(s,1) = mean(results(3).latents(s).acc(data(s).s==1)) - mean(results(3).latents(s).acc(data(s).s==2));
                macc(s,2) = mean(results(3).latents(s).acc(data(s).s==4)) - mean(results(3).latents(s).acc(data(s).s==5));
            end
            
            d = acc(:,1) - acc(:,2);
            err = zeros(2);
            err(:,1) = std(d)./sqrt(length(d));
            m(:,1) = mean(acc);
            
            d = macc(:,1) - macc(:,2);
            err(:,2) = std(d)./sqrt(length(d));
            m(:,2) = mean(acc);
            
            barerrorbar(m',err');
            set(gca,'XTickLabel',{'Data' 'Model'},'FontSize',25,'XLim',[0.5 2.5],'YLim',[0 0.55]);
            ylabel('Go bias','FontSize',25);
            
            [~,p,~,stat] = ttest(acc(:,1),acc(:,2));
            d = mean(acc(:,1)-acc(:,2))./std(acc(:,1)-acc(:,2));
            disp(['go bias (experiment 2): t(',num2str(stat.df),') = ',num2str(stat.tstat),', p = ',num2str(p),', d = ',num2str(d)]);
            
        case 'weight1'
            
            % Pavlovian weight across subjects
            
            data = load_data('data1.csv');
            load results1
            
            for s=1:length(data)
                w(s,:) = [results(1).x(s,2) mean(results(3).latents(s).w)];
                c(s,1) = data(s).cond(1);
            end
            
            plot(w(c==1,1),w(c==1,2),'+b','LineWidth',4,'MarkerSize',10); hold on;
            plot(w(c==2,1),w(c==2,2),'+r','LineWidth',4,'MarkerSize',10);
            legend({'Low control' 'High control'},'FontSize',25,'Location','Northwest');
            set(gca,'XLim',[-0.1 1.1],'YLim',[-0.1 1.1],'XTick',[0 0.5 1],'YTick',[0 0.5 1],'FontSize',25);
            axis square
            ylabel('Average adaptive weight','FontSize',25);
            xlabel('Fixed weight','FontSize',25);
            set(gcf,'Position',[200 200 500 500])
            axis square
            
            [r,p] = corr(w(:,1),w(:,2));
            disp(['r = ',num2str(r),', p = ',num2str(p)]);
            
            data = load_data('data2.csv');
            load results2
            
            clear w
            for s=1:length(data)
                w(s,:) = [results(1).x(s,2) mean(results(3).latents(s).w)];
            end
            
            [r,p] = corr(w(:,1),w(:,2));
            disp(['r = ',num2str(r),', p = ',num2str(p)]);
            
        case 'weight_gobias1'
            
            % Plot go bias as a function of weight quantile in experiment 1
            
            data = load_data('data1.csv');
            load results1
            
            q = linspace(0,1,6);
            for s=1:length(data)
                w = results(3).latents(s).w;
                wq = quantile(w,q);
                for j = 1:length(q)-1
                    ix = w>wq(j) & w<=wq(j+1);
                    go = data(s).s==1;
                    nogo = data(s).s==2;
                    gobias(s,j) = mean(data(s).acc(ix&go)) - mean(data(s).acc(ix&nogo));
                    model_gobias(s,j) = mean(results(3).latents(s).acc(ix&go)) - mean(results(3).latents(s).acc(ix&nogo)); 
                end
            end
            
            q = q(1:end-1) + diff(q)/2;
            [err,m] = wse(gobias);
            errorbar(q,m,err,'ok','LineWidth',4,'MarkerSize',12,'MarkerFaceColor','k');
            hold on; plot(q,nanmedian(model_gobias),'-r','LineWidth',4);
            set(gca,'XLim',[0 1],'FontSize',25,'YLim',[0 0.55]);
            xlabel('Adaptive weight','FontSize',25);
            ylabel('Go bias','FontSize',25);
            legend({'Data' 'Model'}','FontSize',25,'Location','NorthWest');
            
            [~,p,~,stat] = ttest(gobias(:,1),gobias(:,end));
            d = nanmean(gobias(:,1)-gobias(:,end))./nanstd(gobias(:,1)-gobias(:,end));
            disp(['weight (experiment 1): t(',num2str(stat.df),') = ',num2str(stat.tstat),', p = ',num2str(p),', d = ',num2str(d)]);
            
        case 'weight_gobias2'
            
            % Plot go bias as a function of weight quantile in experiment 2
            
            data = load_data('data2.csv');
            if nargin < 2; load results2; end
            
            q = linspace(0,1,6);
            for s=1:length(data)
                w = results(3).latents(s).w;
                wq = quantile(w,q);
                for j = 1:length(q)-1
                    ix = w>wq(j) & w<=wq(j+1);
                    go = data(s).s==1;
                    nogo = data(s).s==2;
                    gobias(s,j) = mean(data(s).acc(ix&go)) - mean(data(s).acc(ix&nogo));
                    model_gobias(s,j) = mean(results(3).latents(s).acc(ix&go)) - mean(results(3).latents(s).acc(ix&nogo)); 
                end
            end
            
            q = q(1:end-1) + diff(q)/2;
            [err,m] = wse(gobias);
            errorbar(q,m,err,'ok','LineWidth',4,'MarkerSize',12,'MarkerFaceColor','k');
            hold on; plot(q,nanmedian(model_gobias),'-r','LineWidth',4);
            set(gca,'XLim',[0 1],'FontSize',25,'YLim',[0 0.55]);
            xlabel('Adaptive weight','FontSize',25);
            ylabel('Go bias','FontSize',25);
            
            [~,p,~,stat] = ttest(gobias(:,1),gobias(:,end));
            d = nanmean(gobias(:,1)-gobias(:,end))./nanstd(gobias(:,1)-gobias(:,end));
            disp(['weight (experiment 2): t(',num2str(stat.df),') = ',num2str(stat.tstat),', p = ',num2str(p),', d = ',num2str(d)]);
            
        case 'simulation'
            
            %rng(3);
            data = load_data('data1.csv');
            load results1;
            param = median(results(3).x);
            data = data(randperm(length(data)));
            simdata(1) = sim_adaptive(param,data(1).s,[0.25 0.75; 0.75 0.25; 0.5 0.5]);
            simdata(2) = sim_adaptive(param,data(2).s,[0.25 0.75; 0.75 0.25; 0.2 0.8]);
            
            w = [simdata(1).w simdata(2).w];
            
            subplot(1,2,1);
            plot(w,'LineWidth',4);
            legend({'Low control' 'High control'},'FontSize',25,'Location','Best');
            set(gca,'XLim',[0 121],'FontSize',25,'YLim',[-0.1 1.1]);
            ylabel('Adaptive weight','FontSize',25);
            xlabel('Trial','FontSize',25);
            
            subplot(1,2,2);
            go_bias(1) = mean(simdata(1).acc(simdata(1).s==1)) - mean(simdata(1).acc(simdata(1).s==2));
            go_bias(2) = mean(simdata(2).acc(simdata(2).s==1)) - mean(simdata(2).acc(simdata(2).s==2));
            bar(go_bias);
            set(gca,'XTickLabel',{'Low control' 'High control'},'FontSize',25,'XLim',[0.5 2.5],'YLim',[0 0.4]);
            ylabel('Go bias','FontSize',25);
            set(gcf,'Position',[200 200 900 400])
            
        case 'bias_variance'
            
            data = load_data('data1.csv');
            load results1
            b = cell(1,2); v = cell(1,2);
            b2 = cell(1,2); v2 = cell(1,2);
            for s = 1:length(data)
                c = data(s).cond(1);
                ix = data(s).s<3;
                go = double(data(s).s(ix)==1);
                a = double(data(s).a(ix)==2);
                b{c}(end+1) = mean(a - go);
                v{c}(end+1) = mean((a-mean(a)).^2);
                a = 1 - results(3).latents(s).P(ix);
                b2{c}(end+1) = mean(a - go);
                v2{c}(end+1) = mean((a-mean(a)).^2);
            end
            
            [~,p,~,stat] = ttest2(b{1},b{2});
            d = (mean(b{1})-mean(b{2}))./sqrt(0.5*(var(b{1})+var(b{2})));
            disp(['Bias, experiment 1: t(',num2str(stat.df),') = ',num2str(stat.tstat),', p = ',num2str(p),', d = ',num2str(d)]);
            [~,p,~,stat] = ttest2(v{1},v{2});
            d = (mean(v{1})-mean(v{2}))./sqrt(0.5*(var(v{1})+var(v{2})));
            disp(['Variance, experiment 1: t(',num2str(stat.df),') = ',num2str(stat.tstat),', p = ',num2str(p),', d = ',num2str(d)]);
            
            b_m = [mean(b{1}) mean(b{2})];
            b_err = [std(b{1})./sqrt(length(b{1})) std(b{2})./sqrt(length(b{2}))];
            v_err = [std(v{1})./sqrt(length(v{1})) std(v{2})./sqrt(length(v{2}))];
            v_m = [mean(v{1}) mean(v{2})];
            
            subplot(2,2,1);
            barerrorbar(b_m',b_err');
            mytitle('a: Experiment 1','Left','FontSize',25,'FontWeight','Bold');
            set(gca,'XTickLabel',{'Low control' 'High control'},'FontSize',25,'XLim',[0.5 2.5],'YLim',[0 0.2]);
            ylabel('Bias','FontSize',25);
            subplot(2,2,2);
            barerrorbar(v_m',v_err');
            set(gca,'XTickLabel',{'Low control' 'High control'},'FontSize',25,'XLim',[0.5 2.5],'YLim',[0.2 0.25]);
            ylabel('Variance','FontSize',25);
            
            data = load_data('data2.csv');
            load results2
            b = []; v = [];
            b2 = []; v2 = [];
            for s = 1:length(data)
                for c = 1:2
                    if c==1
                        ix = data(s).s<3;
                        go = double(data(s).s(ix)==1);
                    else
                        ix = data(s).s<6 & data(s).s>3;
                        go = double(data(s).s(ix)==4);
                    end
                    a = double(data(s).a(ix)==2);
                    b(s,c) = nanmean(a - go);
                    v(s,c) = nanmean((a-mean(a)).^2);
                    a = 1 - results(3).latents(s).P(ix);
                    b2(s,c) = mean(a - go);
                    v2(s,c) = mean((a-mean(a)).^2);
                end
            end
            
            [~,p,~,stat] = ttest(b(:,1),b(:,2));
            d = mean(b(:,1)-b(:,2))./std(b(:,1)-b(:,2));
            disp(['Bias, experiment 2: t(',num2str(stat.df),') = ',num2str(stat.tstat),', p = ',num2str(p),', d = ',num2str(d)]);
            [~,p,~,stat] = ttest(v(:,1),v(:,2));
            d = mean(v(:,1)-v(:,2))./std(v(:,1)-v(:,2));
            disp(['Variance, experiment 2: t(',num2str(stat.df),') = ',num2str(stat.tstat),', p = ',num2str(p),', d = ',num2str(d)]);
            
            [b_err,b_m] = wse(b);
            [v_err,v_m] = wse(v);

            subplot(2,2,3);
            barerrorbar(b_m',b_err');
            mytitle('b: Experiment 2','Left','FontSize',25,'FontWeight','Bold');
            set(gca,'XTickLabel',{'Low control' 'High control'},'FontSize',25,'XLim',[0.5 2.5],'YLim',[0 0.2]);
            ylabel('Bias','FontSize',25);
            subplot(2,2,4);
            barerrorbar(v_m',v_err');
            set(gca,'XTickLabel',{'Low control' 'High control'},'FontSize',25,'XLim',[0.5 2.5],'YLim',[0.2 0.25]);
            ylabel('Variance','FontSize',25);
            
            set(gcf,'Position',[200 200 900 650]);
            
        case 'sim_bias_variance'
            
            load simdata
            data = load_data('data1.csv');
            
            b_m = []; b_err = []; v_m = []; v_err = [];
            for i = 2:3
                b = cell(1,2); v = cell(1,2);
                for s = 1:size(simdata,1)
                    c = data(s).cond(1);
                    ix = simdata(s,i).s<3;
                    go = double(simdata(s,i).s(ix)==1);
                    a = double(simdata(s,i).a(ix)==2);
                    b{c}(end+1) = mean(a - go);
                    v{c}(end+1) = mean((a-mean(a)).^2);
                end
                
                b_m = [b_m; mean(b{1}) mean(b{2})];
                b_err = [b_err; std(b{1})./sqrt(length(b{1})) std(b{2})./sqrt(length(b{2}))];
                v_err = [v_err; std(v{1})./sqrt(length(v{1})) std(v{2})./sqrt(length(v{2}))];
                v_m = [v_m; mean(v{1}) mean(v{2})];
            end
            
            subplot(1,2,1);
            barerrorbar(b_m',b_err');
            set(gca,'XTickLabel',{'Low control' 'High control'},'FontSize',25,'XLim',[0.5 2.5],'YLim',[-0.02 0.5]);
            ylabel('Bias','FontSize',25);
            legend({'Instrumental' 'Pavlovian'},'FontSize',25);
            subplot(1,2,2);
            barerrorbar(v_m',v_err');
            set(gca,'XTickLabel',{'Low control' 'High control'},'FontSize',25,'XLim',[0.5 2.5],'YLim',[0 0.3]);
            ylabel('Variance','FontSize',25);
            set(gcf,'Position',[200 200 1000 450]);
            
        case 'complexity_schematic'
            
            subplot(1,2,1);
            x = linspace(-10,10,1000);
            p = @(x,sd,w,q) w*normpdf(x,0,sd) + (1-w)*unifpdf(x,-q,q);
            plot(x,p(x,10,0.5,1),'-r','LineWidth',4);
            hold on;
            plot(x,p(x,10,0.5,5),'-g','LineWidth',4);
            set(gca,'FontSize',25,'XTick',[],'YTick',[]);
            ylabel('P(model|data)','FontSize',25);
            xlabel('All possible data sets','FontSize',25);
            legend({'Simple model','Complex model'},'FontSize',25,'Location','NorthOutside','Box','off');
            
            subplot(1,2,2);
            x = linspace(0,10,1000);
            y = exp(x);
            plot(x,y,'-k','LineWidth',4);
            hold on;
            plot(x,fliplr(y),'-','LineWidth',4,'Color',[0.5 0.5 0.5]);
            set(gca,'FontSize',25,'XTick',[],'YTick',[]);
            ylabel('Error','FontSize',25);
            xlabel('Model complexity','FontSize',25);
            legend({'Variance' 'Bias^2'},'FontSize',25,'Location','NorthOutside','Box','off')
            
            set(gcf,'Position',[200 200 1000 400])
            
        case 'learning_curves'
            
            data = load_data('data1.csv');
            L = {'Exp 1: Go-to-Win' 'Exp 1: No-Go-To-Win' 'Exp 1: Decoy'; 'Exp 2: Go-to-Win' 'Exp 2: No-Go-To-Win' 'Exp 2: Decoy'};
            
            pGo = cell(1,2); n = [0 0];
            for s = 1:length(data)
                c = data(s).cond(1);
                n(c) = n(c) + 1;
                for i = 1:3
                    pGo{c}(n(c),:,i) = smooth(data(s).a(data(s).s==i)==2);
                end
            end
            
            for i = 1:3
                subplot(2,3,i);
                for c = 1:2
                    m(:,c) = squeeze(nanmean(pGo{c}(:,:,i)));
                    err(:,c) = squeeze(nanstd(pGo{c}(:,:,i))./sqrt(n(c)));
                end
                myeb(m,err);
                set(gca,'FontSize',20,'YLim',[0 1]);
                ylabel('P(Go)','FontSize',25);
                xlabel('Trial','FontSize',25);
                title(L{1,i},'FontSize',25,'FontWeight','Bold');
                if i==1
                    legend({'Low control' 'High control'},'FontSize',20,'Box','off','Location','Best');
                end
            end
            
            data = load_data('data2.csv');
            
            pGo = cell(1,2);
            for s = 1:length(data)
                for i = 1:3
                    pGo{1}(s,:,i) = smooth(data(s).a(data(s).s==i)==2);
                    pGo{2}(s,:,i) = smooth(data(s).a(data(s).s==(i+3))==2);
                end
            end
            
            for i = 1:3
                subplot(2,3,i+3);
                for c = 1:2
                    n(c) = size(pGo{c},1);
                    m(:,c) = squeeze(nanmean(pGo{c}(:,:,i)));
                    err(:,c) = squeeze(nanstd(pGo{c}(:,:,i))./sqrt(n(c)));
                    set(gca,'FontSize',20,'YLim',[0 1]);
                    ylabel('P(Go)','FontSize',25);
                    xlabel('Trial','FontSize',25);
                    title(L{2,i},'FontSize',25,'FontWeight','Bold');
                end
                myeb(m,err);
                
                if i==1
                    legend({'Low control' 'High control'},'FontSize',20,'Box','off','Location','Best');
                end
            end
            
            set(gcf,'Position',[200 200 1200 800])
            
        case 'gobias_dynamic'
            
            data = load_data('data1.csv');
            load results1
            
            gobias = cell(2,2);
            for s = 1:length(data)
                c = data(s).cond(1);
                go = data(s).acc(data(s).s==1);
                nogo = data(s).acc(data(s).s==2);
                gobias{c,1}(end+1,:) = smooth(go-nogo);
                b(s,:) = glmfit(log(1:length(go))',go' - nogo');
                go = results(3).latents(s).acc(data(s).s==1);
                nogo = results(3).latents(s).acc(data(s).s==2);
                gobias{c,2}(end+1,:) = smooth(go-nogo);
            end
            
            [~,p,~,stat] = ttest(b(:,2));
            d = mean(b(:,2))./std(b(:,2));
            disp(['Log trial effect (Experiment 1): t(',num2str(stat.df),') = ',num2str(stat.tstat),', p = ',num2str(p),', d = ',num2str(d)])
            
            for j = 1:2
                for c = 1:2
                    [err(:,c),m(:,c)] = wse(gobias{c,j});
                end
                subplot(2,2,j);
                myeb(m,err);
                set(gca,'FontSize',20,'YLim',[0 1]);
                ylabel('Go bias','FontSize',25);
                xlabel('Trial','FontSize',25);
                if j==1; legend({'Low control' 'High control'},'FontSize',25); end
                
                if j ==1
                    title('Experiment 1: data','FontSize',25,'FontWeight','Bold');
                else
                    title('Experiment 1: model','FontSize',25,'FontWeight','Bold');
                end
            end
            
            data = load_data('data2.csv');
            load results2
            
            gobias = cell(2,2);
            clear b
            for s = 1:length(data)
                go = data(s).acc(data(s).s==1);
                nogo = data(s).acc(data(s).s==2);
                gobias{1,1}(end+1,:) = smooth(go-nogo);
                go = data(s).acc(data(s).s==4);
                nogo = data(s).acc(data(s).s==5);
                gobias{2,1}(end+1,:) = smooth(go-nogo);
                b(s,:) = glmfit(log(1:length(go))',go' - nogo');
                
                go = results(3).latents(s).acc(data(s).s==1);
                nogo = results(3).latents(s).acc(data(s).s==2);
                gobias{1,2}(end+1,:) = smooth(go-nogo);
                go = results(3).latents(s).acc(data(s).s==4);
                nogo = results(3).latents(s).acc(data(s).s==5);
                gobias{2,2}(end+1,:) = smooth(go-nogo);
            end
            
            [~,p,~,stat] = ttest(b(:,2));
            d = mean(b(:,2))./std(b(:,2));
            disp(['Log trial effect (Experiment 2): t(',num2str(stat.df),') = ',num2str(stat.tstat),', p = ',num2str(p),', d = ',num2str(d)])
            
            clear err m
            for j = 1:2
                for c = 1:2
                    [err(:,c),m(:,c)] = wse(gobias{c,j});
                end
                subplot(2,2,j+2);
                myeb(m,err);
                set(gca,'FontSize',20,'YLim',[0 1]);
                ylabel('Go bias','FontSize',25);
                xlabel('Trial','FontSize',25);
                if j ==1
                    title('Experiment 2: data','FontSize',25,'FontWeight','Bold');
                else
                    title('Experiment 2: model','FontSize',25,'FontWeight','Bold');
                end
            end
            
            data = load_data('data2.csv');
            
            set(gcf,'Position',[200 200 1000 800])
            
        case 'instrumental_bias'
            
            load results1
            x = results(3).x(:,2);
            [~,p,~,stat] = ttest(x);
            d = mean(x-0.5)./std(x);
            disp(['Instrumental bias (Experiment 1): t(',num2str(stat.df),') = ',num2str(stat.tstat),', p = ',num2str(p),', d = ',num2str(d)])
            
            load results2
            x = results(3).x(:,2);
            [~,p,~,stat] = ttest(x);
            d = mean(x-0.5)./std(x);
            disp(['Instrumental bias (Experiment 2): t(',num2str(stat.df),') = ',num2str(stat.tstat),', p = ',num2str(p),', d = ',num2str(d)])
    end
    
end

function H = myeb(Y,varargin)
    %
    % myeb(Y,varargin);
    %
    % This function makes nice coloured, shaded error bars. Exactly what
    % it does depends on Y, and on whether you give it one or two inputs.
    %
    % If you only pass it Y, and no other arguments, it assuemd you're
    % giving it raw data.
    %
    %		myeb(Raw_Data)
    %
    % 	.) if Y is 2D array, it will then plot mean(Y) with errorbars given
    % 	by std(Y). In this case there is only one mean vector with its
    % 	errorbars.
    %
    %	.) if Y is 3D array, it will plot size(Y,3) lines with the
    %	associated errorbars. Line k will be mean(Y(:,:,k)) with errorbars
    %	given by std(Y(:,:,k))
    %
    % If you pass it 2 arguments, each has to be at most 2D.
    %
    %		myeb(mu,std)
    %
    % 	.) if mu and std are 1D, it just plots one line given by mu with a
    % 	shaded region given by std.
    %
    %	.) if mu and std are 2D, it will plot size(Y,2) lines in the
    %	standard sequence of colours; each line mu(:,k) will have a shaded
    %	region in the same colour, but less saturated given by std(:,k)
    
    
    %col=[0 0 1; 0 .5 0; 1 0 0; 0 1 1; 1 0 1; 1 .5 0; 1 .5 1];
    %col=[0.8 0.5 0; 0 0 1; 0 .5 0; 1 0 0; 1 0 1; 1 .5 0; 1 .5 1];
    %col = linspecer(size(Y,2));
    %col = col([1 end],:);
    col = linspecer(2);
    ccol=col+.5; ccol(ccol>1)=1;
    
    if isempty(varargin)
        
        if length(size(Y))==2
            m=mean(Y);
            s=std(Y);
            ind1=1:length(m);
            ind2=ind1(end:-1:1);
            %hold on; h=fill([ind1 ind2],[m-s m(ind2)+s(ind2)],.6*ones(1,3));
            hold on; h=fill([ind1 ind2],[m-s m(ind2)+s(ind2)],.6*[1 0 0]);
            set(h,'edgecolor',.6*[1 0 0])
            H = plot(ind1,m,'LineWidth',4);
            hold off
        elseif length(size(Y))>2
            cla; hold on;
            ind1=1:size(Y,2);
            ind2=ind1(end:-1:1);
            if size(Y,3)>8; col=jet(size(Y,3));ccol=col+.8; ccol(ccol>1)=1;end
            for k=1:size(Y,3)
                m=mean(Y(:,:,k));
                s=std(Y(:,:,k));
                h=fill([ind1 ind2],[m-s m(ind2)+s(ind2)],ccol(k,:));
                set(h,'edgecolor',ccol(k,:))
            end
            for k=1:size(Y,3)
                m=mean(Y(:,:,k));
                s=std(Y(:,:,k));
                H = plot(ind1,m,'LineWidth',4,'color',col(k,:));
            end
            hold off
        end
        
    elseif length(varargin)==1
        
        m=Y;
        s=varargin{1};
        if length(size(Y))>2; error;
        elseif min(size(Y))==1
            if size(m,1)>1; m=m';s=s';end
            ind1=1:length(m);
            ind2=ind1(end:-1:1);
            hold on; h=fill([ind1 ind2],[m-s m(ind2)+s(ind2)],.6*ones(1,3));
            set(h,'edgecolor',.6*ones(1,3))
            H = plot(ind1,m,'LineWidth',4);
            hold off
        else
            
            ind1=(1:size(Y,1));
            ind2=ind1(end:-1:1);
            cla; hold on;
            for k=1:size(Y,2)
                H = plot(ind1,m(:,k)','LineWidth',4,'color',col(k,:));
            end
            if size(Y,2)>8; col=jet(size(Y,2));ccol=col+.8; ccol(ccol>1)=1;end
            for k=1:size(Y,2)
                mm=m(:,k)';
                ss=s(:,k)';
                h=fill([ind1 ind2],[mm-ss mm(ind2)+ss(ind2)],ccol(k,:));
                set(h,'edgecolor',ccol(k,:));
            end
            for k=1:size(Y,2)
                mm=m(:,k)';
                H = plot(ind1,mm,'LineWidth',4,'color',col(k,:));
            end
            hold off
        end
        
    elseif length(varargin)==2
        
        m=Y;
        s=varargin{1};
        ix = varargin{2};
        if length(size(Y))>2
            error;
        elseif min(size(Y))==1
            if size(m,1)>1; m=m';s=s';end
            ind1=1:length(m);
            ind2=ind1(end:-1:1);
            hold on; h=fill([ix(ind1) ix(ind2)],[m-s m(ind2)+s(ind2)],.6*ones(1,3));
            set(h,'edgecolor',.6*ones(1,3))
            H = plot(ix(ind1),m,'LineWidth',4);
            hold off
        else
            ind1=(1:size(Y,1));
            ind2=ind1(end:-1:1);
            cla; hold on;
            if size(Y,2)>8; col=jet(size(Y,2));ccol=col+.8; ccol(ccol>1)=1;end
            for k=1:size(Y,2)
                mm=m(:,k)';
                ss=s(:,k)';
                h=fill([ix(ind1) ix(ind2)],[mm-ss mm(ind2)+ss(ind2)],ccol(k,:));
                set(h,'edgecolor',ccol(k,:))
            end
            for k=1:size(Y,2)
                mm=m(:,k)';
                ss=s(:,k)';
                H = plot(ix(ind1),mm,'LineWidth',4,'color',col(k,:));
            end
            hold off
        end
    end
end

function [se, m] = wse(X,dim)
    
    % Within-subject error, following method of Cousineau (2005).
    %
    % USAGE: [se, m] = wse(X,dim)
    %
    % INPUTS:
    %   X - [N x D] data with N observations and D subjects
    %   dim (optional) - dimension along which to compute within-subject
    %   variance (default: 2)
    %
    % OUTPUTS:
    %   se - [1 x D] within-subject standard errors
    %   m - [1 x D] means
    %
    % Sam Gershman, June 2015
    
    if nargin < 2; dim = 2; end
    m = squeeze(nanmean(X));
    X = bsxfun(@minus,X,nanmean(X,dim));
    N = sum(~isnan(X));
    se = bsxfun(@rdivide,nanstd(X),sqrt(N));
    se = squeeze(se);
end