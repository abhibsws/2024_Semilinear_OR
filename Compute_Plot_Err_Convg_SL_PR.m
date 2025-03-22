eqn='SL_PR'; TC=1; tf=1.2;

Lambda = -10.^(1:7);
NT = round(logspace(0.5, 3, 16));  
% DIRK-(s,p,q) scheme
S=[5,5,8,7,10]; 
P=[4,5,4,4,5]; 
Q=[1,1,3,4,4]; 
SchNo=[1,2,3,4,5];


% Compute
DTS = cell(1,length(Lambda)); U_ERR = cell(1,length(Lambda));
tic
for k = 1:length(Lambda)
    lambda = Lambda(k);
    for i = 1:length(S)
        s=S(i); p=P(i); q=Q(i); scheme_no=SchNo(i); 
        for j = 1:length(NT)
            dt = (1/NT(j)); nt = ceil(tf/dt);dt = tf/nt;
            uerr = DIRK_SL_PR_Eqn(TC,tf,s,p,q,scheme_no,dt,nt,lambda);
            DTS{1,k}(i,j) = dt; U_ERR{1,k}(i,j) = uerr;
        end
    end
end
toc

% Create a folder to save figures
save_fig = 1;
if save_fig
    foldername_fig = sprintf('Figures');
    if exist(foldername_fig,'dir')==0,mkdir(foldername_fig);end
end

%
dts = DTS{1,1}(1,:);
% Plot
Cref = {[0.5,0.5,0.5]}; ms = 12; fs = 25;
numColors = length(Lambda);
% Define different shades of blue
% colors = [linspace(0, 0, numColors)', linspace(0.2, 0.8, numColors)', linspace(0.5, 1, numColors)']; 
colors = winter(length(Lambda));
% Define at least 8 different markers
markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '*', '+'}; linS = {':','--','-'};

yu = [5e-4;5e-3;1e-5;5e-6;3e-7];
yl = [3e-14;5e-16;1e-16;1e-16;1e-16];

or = [1,2,4;1,2,5;2,3,4;2,3,4;2,4,5];

st = [1,1,1;1,1,1;1,1,6;1,1,6;1,1,6]; 
en = [8,16,8;6,16,10;6,6,12;6,8,14;6,6,12]; 

coeff = [4e-8,5e-3,3e-3;
         2e-2,5e-2,1e-2;
         6e-5,5e-8,1e-2;
         5e-5,2e-8,1e-3;
         3e-6,3e-9,1e-2]; 

% Reference methods
for i = 1:length(S)
    %---------for minimum white space
    % Ensure figure remains valid
    fig = figure(i);
    % Set figure size
    set(fig, 'Units', 'pixels', 'Position', [50, 50, 600, 600]);
    % Adjust axis to remove extra margins
    ax = gca;
    set(ax, 'LooseInset', [0, 0, 0, 0]);
    %---------for minimum white space
    legendEntries = {}; s = S(i); p = P(i); q = Q(i);
    for k = 1:length(Lambda)
        lambda = Lambda(k);
        
        loglog(DTS{1,k}(i,:), U_ERR{1,k}(i,:), ...
            ['-', markers{k}], 'Color', colors(k, :), 'MarkerSize', ms, 'LineWidth', 3); 
        legendEntries{end+1} = sprintf('$\\lambda = -10^{%d}$', log10(abs(lambda)));
        hold on         
    end
    % reference slopes
    for j =1:3
        loglog(dts(st(i,j):en(i,j)),coeff(i,j)*dts(st(i,j):en(i,j)).^or(i,j),'LineStyle', linS{j},'color',Cref{1})
        legendEntries{end+1} = sprintf('Slope %d',or(i,j));
        hold on
    end
    
    xlim([dts(1,end),dts(1,1)])
    ylim([yl(i),yu(i)])
    if q == 1
        title(sprintf('$(%d,%d,%d)$', s, p, q), 'Interpreter', 'latex'); 
        xlabel('h', 'Position', get(gca, 'XLabel').Position + [0, 4e-17, 0])
        ylabel('Error', 'Position', get(gca, 'YLabel').Position + [-2e-5, 1e-7, 0])
    else
        title(sprintf('$(\\mathbf{%d},\\mathbf{%d},\\mathbf{%d})$', s, p, q), 'Interpreter', 'latex'); 
        xlabel('h', 'Position', get(gca, 'XLabel').Position + [0, 4e-17, 0])
        ylabel('Error', 'Position', get(gca, 'YLabel').Position + [-2e-5, 1e-12, 0])
    end
    grid minor
    set(gca,'FontSize',fs)
    legend(legendEntries, 'Interpreter', 'latex','NumColumns',2,'Location','best','Box','off','FontSize',fs)
    % Save as pdf
    figure_name = sprintf('Figures/%s_TC%d_tf%1.1f_Convg_s%dp%dq%d.pdf',eqn,TC,tf,s,p,q);
    print(gcf,figure_name,'-dpdf','-r100','-bestfit')
end


%%=============================Old code==================================%%
% eqn='SL_PR'; TC=1; tf=1.2;
% 
% Lambda = [-1,-10,-10^2,-10^3,-10^4];
% % Lambda = [-10^1,-10^4,-10^8];
% NT = 2.^(0:.5:9);  
% % DIRK-(s,p,q) scheme
% S=[8,5,7,5,10,5]; 
% P=[4,4,4,4,5,5]; 
% Q=[3,1,4,1,4,1]; 
% SchNo=[3,1,4,1,5,2];
% 
% 
% % Compute
% DTS = cell(1,length(Lambda)); U_ERR = cell(1,length(Lambda));
% tic
% for k = 1:length(Lambda)
%     lambda = Lambda(k);
%     for i = 1:length(S)
%         s=S(i); p=P(i); q=Q(i); scheme_no=SchNo(i); 
%         for j = 1:length(NT)
%             dt = (1/NT(j)); nt = ceil(tf/dt);dt = tf/nt;
%             uerr = DIRK_SL_PR_Eqn(TC,tf,s,p,q,scheme_no,dt,nt,lambda);
%             DTS{1,k}(i,j) = dt; U_ERR{1,k}(i,j) = uerr;
%         end
%     end
% end
% toc
% 
% % Create a folder to save figures
% save_fig = 1;
% if save_fig
%     foldername_fig = sprintf('Figures');
%     if exist(foldername_fig,'dir')==0,mkdir(foldername_fig);end
% end
% 
% %
% dts = DTS{1,1}(1,:);
% 
% % Plot
% legendEntries = {}; C = {'b','r','m','k','g','y'}; Cref = {[0.5,0.5,0.5]};
% linS = {'-','--',':'}; Mar = {'o','s','+','*','d'}; ms = 10; fs = 26;
% figure(1)
% set(gcf,'position',[0 0 800 1000])
% set(gca, 'LooseInset', max(get(gca, 'TightInset'), 0)); % remove white space
% set(0,'DefaultLineLineWidth',3);
% for k = 1:length(Lambda)
%     lambda = Lambda(k);
%     for i = 1:2
%         s=S(i); p=P(i); q=Q(i); scheme_no=SchNo(i); 
%         if mod(i, 2) == 1
%             loglog(DTS{1,k}(i,:),U_ERR{1,k}(i,:),'s-','color',C{k},'marker',Mar{k},'MarkerSize',ms);
%             hold on
%             legendEntries{end+1} = sprintf('$\\lambda = -10^{%d}: (\\mathbf{%d},\\mathbf{%d},\\mathbf{%d})$', log10(abs(lambda)), s, p, q);            
%         else 
%             loglog(DTS{1,k}(i,:),U_ERR{1,k}(i,:),'--','color',C{k},'marker',Mar{k},'MarkerSize',ms);
%             hold on
%             legendEntries{end+1} = sprintf('$\\lambda = -10^{%d}: (%d,%d,%d)$',log10(abs(lambda)), s, p, q);
%         end
%     end
% end
% 
% % reference slope
% sl1 = [2,3,4];
% st11 = 1; en11 = 5; Coef11 = 5e-3; 
% st12 = 1; en12 = 5; Coef12 = 5e-7;
% st13 = 7; en13 = 10; Coef13 = 8e-3;
% 
% loglog(dts(st11:en11),Coef11*dts(st11:en11).^sl1(1),':','color',Cref{1})
% hold on
% legendEntries{end+1} = sprintf('Slope %d',sl1(1));
% loglog(dts(st12:en12),Coef12*dts(st12:en12).^sl1(2),'--','color',Cref{1})
% hold on
% legendEntries{end+1} = sprintf('Slope %d',sl1(2));
% loglog(dts(st13:en13),Coef13*dts(st13:en13).^sl1(3),'-','color',Cref{1})
% hold on
% legendEntries{end+1} = sprintf('Slope %d',sl1(3));
% 
% xlim([dts(1,end),dts(1,1)])
% ylim([5e-17,5e-3])
% xlabel('\Delta t');
% ylabel('Error', 'Position', get(gca, 'YLabel').Position + [-1e-5, 1e-7, 0])
% grid minor
% set(gca,'FontSize',fs)
% 
% % Set legend
% legend(legendEntries, 'Interpreter', 'latex','NumColumns',2,'Location','northoutside','Box','off','FontSize',fs)
% 
% % Save as pdf
% figure_name = sprintf('Figures/%s_TC%d_tf%1.1f_Convg_s%dp%dq%d.pdf',eqn,TC,tf,S(1),P(1),Q(1));
% print(gcf,figure_name,'-dpdf','-r100','-bestfit')
% 
% % Plot
% figure(2)
% set(gcf,'position',[0 0 800 1000])
% set(gca, 'LooseInset', max(get(gca, 'TightInset'), 0)); % remove white space
% set(0,'DefaultLineLineWidth',3);
% legendEntries = {}; 
% for k = 1:length(Lambda)
%     lambda = Lambda(k);
%     for i = 3:4
%         s=S(i); p=P(i); q=Q(i); scheme_no=SchNo(i); 
%         if mod(i, 2) == 1
%             loglog(DTS{1,k}(i,:),U_ERR{1,k}(i,:),'s-','color',C{k},'marker',Mar{k},'MarkerSize',ms);
%             hold on
%             legendEntries{end+1} = sprintf('$\\lambda = -10^{%d}: (\\mathbf{%d},\\mathbf{%d},\\mathbf{%d})$', log10(abs(lambda)), s, p, q);        
%         else 
%             loglog(DTS{1,k}(i,:),U_ERR{1,k}(i,:),'--','color',C{k},'marker',Mar{k},'MarkerSize',ms);
%             hold on
%             legendEntries{end+1} = sprintf('$\\lambda = -10^{%d}: (%d,%d,%d)$',log10(abs(lambda)), s, p, q);
%         end
%     end
% end
% 
% % reference slope
% sl2 = [2,3,4];
% st21 = 1; en21 = 6; Coef21 = 5e-3; 
% st22 = 1; en22 = 5; Coef22 = 7e-7;
% st23 = 5; en23 = 10; Coef23 = 2e-3;
% 
% loglog(dts(st21:en21),Coef21*dts(st21:en21).^sl2(1),':','color',Cref{1})
% hold on
% legendEntries{end+1} = sprintf('Slope %d',sl2(1));
% loglog(dts(st22:en22),Coef22*dts(st22:en22).^sl2(2),'--','color',Cref{1})
% hold on
% legendEntries{end+1} = sprintf('Slope %d',sl2(2));
% loglog(dts(st23:en23),Coef23*dts(st23:en23).^sl2(3),'-','color',Cref{1})
% hold on
% legendEntries{end+1} = sprintf('Slope %d',sl2(3));
% 
% xlim([dts(1,end),dts(1,1)])
% ylim([5e-17,5e-3])
% xlabel('\Delta t');
% ylabel('Error', 'Position', get(gca, 'YLabel').Position + [-1e-5, 1e-7, 0])
% grid minor
% set(gca,'FontSize',fs)
% 
% 
% % Set legend
% legend(legendEntries, 'Interpreter', 'latex','NumColumns',2,'Location','northoutside','Box','off','FontSize',fs)
% 
% % Save as pdf
% figure_name = sprintf('Figures/%s_TC%d_tf%1.1f_Convg_s%dp%dq%d.pdf',eqn,TC,tf,S(3),P(3),Q(3));
% print(gcf,figure_name,'-dpdf','-r100','-bestfit')
% 
% % Plot
% figure(3)
% set(gcf,'position',[0 0 800 1000])
% set(gca, 'LooseInset', max(get(gca, 'TightInset'), 0)); % remove white space
% set(0,'DefaultLineLineWidth',3);
% legendEntries = {};
% for k = 1:length(Lambda)
%     lambda = Lambda(k);
%     for i = 5:6
%         s=S(i); p=P(i); q=Q(i); scheme_no=SchNo(i); 
%         if mod(i, 2) == 1
%             loglog(DTS{1,k}(i,:),U_ERR{1,k}(i,:),'s-','color',C{k},'marker',Mar{k},'MarkerSize',ms);
%             hold on
%             legendEntries{end+1} = sprintf('$\\lambda = -10^{%d}: (\\mathbf{%d},\\mathbf{%d},\\mathbf{%d})$', log10(abs(lambda)), s, p, q);          
%         else 
%             loglog(DTS{1,k}(i,:),U_ERR{1,k}(i,:),'--','color',C{k},'marker',Mar{k},'MarkerSize',ms);
%             hold on
%             legendEntries{end+1} = sprintf('$\\lambda = -10^{%d}: (%d,%d,%d)$',log10(abs(lambda)), s, p, q);
%         end
%     end
% end
% 
% % reference slope
% sl3 = [2,4,5];
% st31 = 1; en31 = 6; Coef31 = 1e-2; 
% st32 = 1; en32 = 4; Coef32 = 7e-9;
% st33 = 5; en33 = 9; Coef33 = 6e-3;
% 
% loglog(dts(st31:en31),Coef31*dts(st31:en31).^sl3(1),':','color',Cref{1})
% hold on
% legendEntries{end+1} = sprintf('Slope %d',sl3(1));
% loglog(dts(st32:en32),Coef32*dts(st32:en32).^sl3(2),'--','color',Cref{1})
% hold on
% legendEntries{end+1} = sprintf('Slope %d',sl3(2));
% loglog(dts(st33:en33),Coef33*dts(st33:en33).^sl3(3),'-','color',Cref{1})
% hold on
% legendEntries{end+1} = sprintf('Slope %d',sl3(3));
% 
% xlim([dts(1,end),dts(1,1)])
% ylim([4e-17,8e-2])
% xlabel('\Delta t');
% ylabel('Error', 'Position', get(gca, 'YLabel').Position + [-1e-5, 1e-7, 0])
% grid minor
% set(gca,'FontSize',fs)
% 
% % Set legend
% legend(legendEntries, 'Interpreter', 'latex','NumColumns',2,'Location','northoutside','Box','off','FontSize',fs)
% 
% % Save as pdf
% figure_name = sprintf('Figures/%s_TC%d_tf%1.1f_Convg_s%dp%dq%d.pdf',eqn,TC,tf,S(5),P(5),Q(5));
% print(gcf,figure_name,'-dpdf','-r100','-bestfit')
