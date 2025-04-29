eqn='VanDerPol'; TC=1; tf=0.5;
Epsilon = [1e-3,1e-6];
NT = 2.^(2:1:11);  
% DIRK-(s,p,q) scheme
S=[8,5,10,5]; 
P=[4,4,5,5]; 
Q=[3,1,4,1]; 
SchNo=[3,1,5,2];


% Compute
DTS = cell(1,length(Epsilon)); U_ERR = cell(1,length(Epsilon));
tic
for k = 1:length(Epsilon)
    ep = Epsilon(k);
    for i = 1:length(S)
        s=S(i); p=P(i); q=Q(i); scheme_no=SchNo(i); 
        for j = 1:length(NT)
            dt = (1/NT(j)); nt = ceil(tf/dt);dt = tf/nt;
            uerr = DIRKspqVanDerPol(TC,tf,s,p,q,scheme_no,dt,nt,ep);
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
C = {'b','r','g','k','m'}; Cref = {[0.5,0.5,0.5]};
linS = {'-','--',':'}; Mar = {'o','s','+','*','d'}; ms = 10; fs = 20;

% Plot
%---------for minimum white space
% Ensure figure remains valid
fig = figure(1);
% Set figure size
set(fig, 'Units', 'pixels', 'Position', [50, 50, 500, 500]);
% Adjust axis to remove extra margins
ax = gca;
set(ax, 'LooseInset', [0, 0, 0, 0]);
set(0,'DefaultLineLineWidth',3);
%---------for minimum white space
legendEntries = {};
for k = 1:length(Epsilon)
    exp = -log10(Epsilon(k));
    for i = 1:2
        s=S(i); p=P(i); q=Q(i); scheme_no=SchNo(i); 
        if mod(i, 2) == 1
            loglog(DTS{1,k}(i,:),U_ERR{1,k}(i,:),'s-','color',C{k},'marker',Mar{k},'MarkerSize',ms);
            hold on 
            legendEntries{end+1} = sprintf('$\\epsilon = 10^{-%d}: (\\mathbf{%d},\\mathbf{%d},\\mathbf{%d})$', exp, s, p, q);
        else 
            loglog(DTS{1,k}(i,:),U_ERR{1,k}(i,:),'--','color',C{k},'marker',Mar{k},'MarkerSize',ms);
            hold on
            legendEntries{end+1} = sprintf('$\\epsilon = 10^{-%d}: (%d,%d,%d)$', exp, s, p, q);
        end
    end
end
% reference slope
sl1 = [1,3,4];
st11 = 1; en11 = 6; Coef11 = 1e-3; 
st12 = 4; en12 = 9; Coef12 = 6e-3;
st13 = 3; en13 = 8; Coef13 = 3e-3;

loglog(dts(st11:en11),Coef11*dts(st11:en11).^sl1(1),':','color',Cref{1})
hold on
legendEntries{end+1} = sprintf('Slope %d',sl1(1));
loglog(dts(st12:en12),Coef12*dts(st12:en12).^sl1(2),'--','color',Cref{1})
hold on
legendEntries{end+1} = sprintf('Slope %d',sl1(2));
loglog(dts(st13:en13),Coef13*dts(st13:en13).^sl1(3),'-','color',Cref{1})
hold on
legendEntries{end+1} = sprintf('Slope %d',sl1(3));


xlim([dts(1,end),dts(1,1)])
ylim([8e-15,5e-4])
xlabel('h', 'Position', get(gca, 'XLabel').Position + [-5e-3, 2e-15, 0])
ylabel('Error', 'Position', get(gca, 'YLabel').Position + [0.6e-4, 1e-7, 0])
grid minor
set(gca,'FontSize',fs)
% Tighten axis limits and position
% axis tight
% set(gca, 'Position', [0.1 0.1 0.80 0.80]); % Adjust to reduce margins
% Set legend
legend(legendEntries, 'Interpreter', 'latex','NumColumns',1,'Location','southeast', 'Box', 'off','FontSize',fs)

% Save as pdf
figure_name = sprintf('Figures/%s_TC%d_tf%1.1f_Convg_s%dp%dq%d.pdf',eqn,TC,tf,S(1),P(1),Q(1));
%print(gcf,figure_name,'-dpdf','-r100','-bestfit')
exportgraphics(gcf, figure_name, 'ContentType', 'vector')


% Plot
%---------for minimum white space
% Ensure figure remains valid
fig = figure(2);
% Set figure size
set(fig, 'Units', 'pixels', 'Position', [50, 50, 500, 500]);
% Adjust axis to remove extra margins
ax = gca;
set(ax, 'LooseInset', [0, 0, 0, 0]);
set(0,'DefaultLineLineWidth',3);
%---------for minimum white space
legendEntries = {};
for k = 1:length(Epsilon)
    exp = -log10(Epsilon(k));
    for i = 3:4
        s=S(i); p=P(i); q=Q(i); scheme_no=SchNo(i); 
        if mod(i, 2) == 1
            loglog(DTS{1,k}(i,:),U_ERR{1,k}(i,:),'s-','color',C{k},'marker',Mar{k},'MarkerSize',ms);
            hold on
            legendEntries{end+1} = sprintf('$\\epsilon = 10^{-%d}: (\\mathbf{%d},\\mathbf{%d},\\mathbf{%d})$', exp, s, p, q);
        else 
            loglog(DTS{1,k}(i,:),U_ERR{1,k}(i,:),'--','color',C{k},'marker',Mar{k},'MarkerSize',ms);
            hold on
            legendEntries{end+1} = sprintf('$\\epsilon = 10^{-%d}: (%d,%d,%d)$', exp, s, p, q);
        end
    end
end

% reference slope
sl2 = [2,4,5];
st21 = 1; en21 = 6; Coef21 = 8e-1; 
st22 = 4; en22 = 9; Coef22 = 5e-1;
st23 = 3; en23 = 6; Coef23 = 1e-2;

loglog(dts(st21:en21),Coef21*dts(st21:en21).^sl2(1),':','color',Cref{1})
hold on
legendEntries{end+1} = sprintf('Slope %d',sl2(1));
loglog(dts(st22:en22),Coef22*dts(st22:en22).^sl2(2),'--','color',Cref{1})
hold on
legendEntries{end+1} = sprintf('Slope %d',sl2(2));
loglog(dts(st23:en23),Coef23*dts(st23:en23).^sl2(3),'-','color',Cref{1})
hold on
legendEntries{end+1} = sprintf('Slope %d',sl2(3));

xlim([dts(1,end),dts(1,1)])
ylim([7e-16,2e-0])
xlabel('h', 'Position', get(gca, 'XLabel').Position + [-0.6e-2, 1.5e-16, 0])
ylabel('Error', 'Position', get(gca, 'YLabel').Position + [1e-4, 0.5e-7, 0])
grid minor
set(gca,'FontSize',fs)
% Set legend
legend(legendEntries, 'Interpreter', 'latex','NumColumns',1,'Location','northwest', 'Box', 'off','FontSize',fs)

% Save as pdf
figure_name = sprintf('Figures/%s_TC%d_tf%1.1f_Convg_s%dp%dq%d.pdf',eqn,TC,tf,S(3),P(3),Q(3));
% print(gcf,figure_name,'-dpdf','-r100','-bestfit')
exportgraphics(gcf, figure_name, 'ContentType', 'vector')
