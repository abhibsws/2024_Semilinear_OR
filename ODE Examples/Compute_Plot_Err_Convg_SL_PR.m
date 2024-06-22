eqn='SL_PR'; TC=1; tf=1;

Lambda = [-10,-10^2,-10^3,-10^4,-10^5,-10^6];
NT = 2.^(0:1:10);  
% DIRK-(s,p,q) scheme
S=[6,5,7,5,10,5]; 
P=[4,4,4,4,5,5]; 
Q=[3,1,4,1,4,1]; 
SchNo=[8,2,9,2,10,3];


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
legendEntries = {}; C = {'b','r','g','k','m','y'}; Cref = {[0.5,0.5,0.5]};
linS = {'-','--',':'}; ms = 4; fs = 16;
figure(1)
set(gcf,'position',[0 0 600 800])
set(0,'DefaultLineLineWidth',3);
for k = 1:length(Lambda)
    lambda = Lambda(k);
    for i = 1:2
        s=S(i); p=P(i); q=Q(i); scheme_no=SchNo(i); 
        if mod(i, 2) == 1
            loglog(DTS{1,k}(i,:),U_ERR{1,k}(i,:),'s-','color',C{k},'MarkerSize',ms);
        else 
            loglog(DTS{1,k}(i,:),U_ERR{1,k}(i,:),'--','color',C{k},'MarkerSize',ms);
        end
        hold on
        % Create legend entry
        legendEntries{end+1} = sprintf('$\\lambda = -10^{%d}: (%d,%d,%d)$', k, s, p, q);
    end
end

% reference slope
sl1 = [2,3,4];
st11 = 1; en11 = 6; Coef11 = 1e-3; 
st12 = 5; en12 = 10; Coef12 = 7e-7;
st13 = 5; en13 = 10; Coef13 = 8e-3;

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
ylim([5e-17,1e-2])
xlabel('\Delta t');
ylabel('Error');
grid minor
set(gca,'FontSize',fs)
% Set legend
legend(legendEntries, 'Interpreter', 'latex','NumColumns',3,'Location','northoutside')

% Save as pdf
figure_name = sprintf('Figures/%s_TC%d_tf%1.1f_Convg_s%dp%dq%d.pdf',eqn,TC,tf,S(1),P(1),Q(1));
print(gcf,figure_name,'-dpdf','-r100','-bestfit')

% Plot
figure(2)
set(gcf,'position',[0 0 600 800])
set(0,'DefaultLineLineWidth',3);
legendEntries = {}; 
for k = 1:length(Lambda)
    lambda = Lambda(k);
    for i = 3:4
        s=S(i); p=P(i); q=Q(i); scheme_no=SchNo(i); 
        if mod(i, 2) == 1
            loglog(DTS{1,k}(i,:),U_ERR{1,k}(i,:),'s-','color',C{k},'MarkerSize',ms);
        else 
            loglog(DTS{1,k}(i,:),U_ERR{1,k}(i,:),'--','color',C{k},'MarkerSize',ms);
        end
        hold on
        % Create legend entry
        legendEntries{end+1} = sprintf('$\\lambda = -10^{%d}: (%d,%d,%d)$', k, s, p, q);
    end
end

% reference slope
sl2 = [2,3,4];
st21 = 1; en21 = 6; Coef21 = 1e-3; 
st22 = 5; en22 = 10; Coef22 = 7e-7;
st23 = 5; en23 = 10; Coef23 = 8e-3;

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
ylim([5e-17,1e-2])
xlabel('\Delta t');
ylabel('Error');
grid minor
set(gca,'FontSize',fs)
% Set legend
legend(legendEntries, 'Interpreter', 'latex','NumColumns',3,'Location','northoutside')

% Save as pdf
figure_name = sprintf('Figures/%s_TC%d_tf%1.1f_Convg_s%dp%dq%d.pdf',eqn,TC,tf,S(3),P(3),Q(3));
print(gcf,figure_name,'-dpdf','-r100','-bestfit')

% Plot
figure(3)
set(gcf,'position',[0 0 600 800])
set(0,'DefaultLineLineWidth',3);
legendEntries = {};
for k = 1:length(Lambda)
    lambda = Lambda(k);
    for i = 5:6
        s=S(i); p=P(i); q=Q(i); scheme_no=SchNo(i); 
        if mod(i, 2) == 1
            loglog(DTS{1,k}(i,:),U_ERR{1,k}(i,:),'s-','color',C{k},'MarkerSize',ms);
        else 
            loglog(DTS{1,k}(i,:),U_ERR{1,k}(i,:),'--','color',C{k},'MarkerSize',ms);
        end
        hold on
        % Create legend entry
        legendEntries{end+1} = sprintf('$\\lambda = -10^{%d}: (%d,%d,%d)$', k, s, p, q);
    end
end

% reference slope
sl3 = [2,4,5];
st31 = 1; en31 = 6; Coef31 = 1e-3; 
st32 = 5; en32 = 10; Coef32 = 7e-7;
st33 = 5; en33 = 10; Coef33 = 8e-3;

loglog(dts(st31:en31),Coef31*dts(st31:en31).^sl3(1),':','color',Cref{1})
hold on
legendEntries{end+1} = sprintf('Slope %d',sl3(1));
loglog(dts(st32:en32),Coef22*dts(st32:en32).^sl3(2),'--','color',Cref{1})
hold on
legendEntries{end+1} = sprintf('Slope %d',sl3(2));
loglog(dts(st33:en33),Coef33*dts(st33:en33).^sl3(3),'-','color',Cref{1})
hold on
legendEntries{end+1} = sprintf('Slope %d',sl3(3));

xlim([dts(1,end),dts(1,1)])
ylim([4e-17,1e-1])
xlabel('\Delta t');
ylabel('Error');
grid minor
set(gca,'FontSize',fs)
% Set legend
legend(legendEntries, 'Interpreter', 'latex','NumColumns',3,'Location','northoutside')

% Save as pdf
figure_name = sprintf('Figures/%s_TC%d_tf%1.1f_Convg_s%dp%dq%d.pdf',eqn,TC,tf,S(5),P(5),Q(5));
print(gcf,figure_name,'-dpdf','-r100','-bestfit')
