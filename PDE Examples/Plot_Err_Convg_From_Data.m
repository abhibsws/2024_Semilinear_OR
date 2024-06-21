% Different test problems
EQN = {'SemiLinAdv','VisBurgers','AllenCahn2D'};

% Different initial conditions
IC = [1,1,1];
% Final times
TF = [0.5,1,0.5];
% Spatial order 
sp_or = {5,6,0};
% Number of grid points
M = [1e3,1e3,30];

%plot colors and linestyles
C = {'b','r','g'}; Cref = {[0.5,0.5,0.5]};
linS = {'-','--',':'}; Mar = {'o','s','+'}; ms = 4; fs = 16;

% Create a folder to save figures
save_fig = 1;
if save_fig
    foldername_fig = sprintf('Figures');
    if exist(foldername_fig,'dir')==0,mkdir(foldername_fig);end
end

% loading data
for eq_ind = 1:3
    eqn=EQN{eq_ind};TC=IC(eq_ind); tf=TF(eq_ind); 
    m = M(eq_ind); spatial_order = sp_or{eq_ind};
    % for legends
    switch eqn 
        case 'SemiLinAdv'
            S=[5,7,5,10]; P=[4,4,5,5]; Q=[1,4,1,4]; SchNo=[2,9,3,10];
            filename= sprintf('errdata/%s/ErrConvgData_%s_TC%d_tf%.1f.mat',eqn,eqn,TC,tf);
            load(filename)

            st = [2,2,3,3]; en = [5,5,5,5];
            Cof = [3e-0,2e-0,5e2,1e3]; Sl = [2,4,2,5]; 
            yl = [6e-10,2e-11]; yu = [1e-1,1e-0];
            % Plot
            for i = 1:2
                figure(i)
                set(0,'DefaultLineLineWidth',3);
                loglog(dts(2*i-1,:),U_Err(2*i-1,:),'linestyle',linS{1},'color',C{1},'marker',Mar{1},'MarkerSize',ms)
                hold on
                loglog(dts(2*i,:),U_Err(2*i,:),'linestyle',linS{2},'color',C{2},'marker',Mar{2},'MarkerSize',ms)
                hold on
                loglog(dts(2*i-1,st(2*i-1):en(2*i-1)),Cof(2*i-1)*dts(2*i-1,st(2*i-1):en(2*i-1)).^Sl(2*i-1),'--','color',Cref{1})
                loglog(dts(2*i,st(2*i):en(2*i)),Cof(2*i)*dts(2*i,st(2*i):en(2*i)).^Sl(2*i),'-','color',Cref{1})
                legend('','Location', 'best')
                legend({sprintf('(%d,%d,%d)',S(2*i-1),P(2*i-1),Q(2*i-1)),...
                sprintf('(%d,%d,%d)-SL',S(2*i),P(2*i),Q(2*i)),...
                sprintf('Slope %d',Sl(2*i-1)),...
                sprintf('Slope %d',Sl(2*i))},'NumColumns',2)
                %
                xlim([dts(1,end),dts(1,1)])
                %ylim([yl(i),yu(i)])
                xlabel('\Delta t');
                ylabel('Error');
                grid minor
                set(gca,'FontSize',fs)
                %
                figure_name = sprintf('Figures/%s_TC%d_tf%1.1f_Convg_Order_s%dp%dq%d_m%d_sp_%d.pdf',eqn,TC,tf,S(2*i),P(2*i),Q(2*i),m,spatial_order);
                print(gcf,figure_name,'-dpdf','-r100','-vector')
            end
        case 'VisBurgers'
            S=[5,7,7,5,12,14]; P=[4,4,4,5,5,5]; Q=[1,4,4,1,5,4]; SchNo=[2,4,9,3,6,11];
            filename= sprintf('errdata/%s/ErrConvgData_%s_TC%d_tf%.1f.mat',eqn,eqn,TC,tf);
            load(filename) 

            st = [4,4,4,4,2,2]; en = [6,6,6,6,6,6];
            Cof = [1e-0,1e-1,8e0,3e1,1e-1,1e0]; Sl = [2,3,4,2,4,5]; 
            yl = [5e-11,1e-12]; yu = [1e-2,1e-1]; 
            % Plot
            for i = 1:2
                figure(2+i)
                set(0,'DefaultLineLineWidth',3);
                loglog(dts(3*i-2,:),U_Err(3*i-2,:),'linestyle',linS{1},'color',C{1},'marker',Mar{1},'MarkerSize',ms)
                hold on
                loglog(dts(3*i-1,:),U_Err(3*i-1,:),'linestyle',linS{2},'color',C{2},'marker',Mar{2},'MarkerSize',ms)
                hold on
                loglog(dts(3*i,:),U_Err(3*i,:),'linestyle',linS{3},'color',C{3},'marker',Mar{3},'MarkerSize',ms)
                hold on
                loglog(dts(3*i-2,st(3*i-2):en(3*i-2)),Cof(3*i-2)*dts(3*i-2,st(3*i-2):en(3*i-2)).^Sl(3*i-2),':','color',Cref{1})
                loglog(dts(3*i-1,st(3*i-1):en(3*i-1)),Cof(3*i-1)*dts(3*i-1,st(3*i-1):en(3*i-1)).^Sl(3*i-1),'--','color',Cref{1})
                loglog(dts(3*i,st(3*i):en(3*i)),Cof(3*i)*dts(3*i,st(3*i):en(3*i)).^Sl(3*i),'-','color',Cref{1})
                legend('','Location', 'best')
                legend({sprintf('(%d,%d,%d)',S(3*i-2),P(3*i-2),Q(3*i-2)),...
                sprintf('(%d,%d,%d)',S(3*i-1),P(3*i-1),Q(3*i-1)),...
                sprintf('(%d,%d,%d)-SL',S(3*i),P(3*i),Q(3*i)),...
                sprintf('Slope %d',Sl(3*i-2)), sprintf('Slope %d',Sl(3*i-1)),...
                sprintf('Slope %d',Sl(3*i))},'NumColumns',2)
                %
                xlim([dts(1,end),dts(1,1)])
                %ylim([yl(i),yu(i)])
                xlabel('\Delta t');
                ylabel('Error');
                grid minor
                set(gca,'FontSize',fs)
                %
                figure_name = sprintf('Figures/%s_TC%d_tf%1.1f_Convg_Order_s%dp%dq%d_m%d_sp_%d.pdf',eqn,TC,tf,S(2*i),P(2*i),Q(2*i),m,spatial_order);
                print(gcf,figure_name,'-dpdf','-r100','-vector')
            end
        case 'AllenCahn2D'
            S=[5,7,5,10]; P=[4,4,5,5]; Q=[1,4,1,4]; SchNo=[2,9,3,10];
            filename= sprintf('errdata/%s/ErrConvgData_%s_TC%d_tf%.1f.mat',eqn,eqn,TC,tf);
            load(filename)

            st = [4,4,3,3]; en = [6,6,6,6];
            Cof = [3e-0,2e-0,5e2,1e3]; Sl = [2,4,2,5]; 
            yl = [6e-10,2e-11]; yu = [1e-1,1e-0];
            % Plot
            for i = 1:2
                figure(4+i)
                set(0,'DefaultLineLineWidth',3);
                loglog(dts(2*i-1,:),U_Err(2*i-1,:),'linestyle',linS{1},'color',C{1},'marker',Mar{1},'MarkerSize',ms)
                hold on
                loglog(dts(2*i,:),U_Err(2*i,:),'linestyle',linS{2},'color',C{2},'marker',Mar{2},'MarkerSize',ms)
                hold on
                loglog(dts(2*i-1,st(2*i-1):en(2*i-1)),Cof(2*i-1)*dts(2*i-1,st(2*i-1):en(2*i-1)).^Sl(2*i-1),'--','color',Cref{1})
                loglog(dts(2*i,st(2*i):en(2*i)),Cof(2*i)*dts(2*i,st(2*i):en(2*i)).^Sl(2*i),'-','color',Cref{1})
                legend('','Location', 'best')
                legend({sprintf('(%d,%d,%d)',S(2*i-1),P(2*i-1),Q(2*i-1)),...
                sprintf('(%d,%d,%d)-SL',S(2*i),P(2*i),Q(2*i)),...
                sprintf('Slope %d',Sl(2*i-1)),...
                sprintf('Slope %d',Sl(2*i))},'NumColumns',2)
                %
                xlim([dts(1,end),dts(1,1)])
                %ylim([yl(i),yu(i)])
                xlabel('\Delta t');
                ylabel('Error');
                grid minor
                set(gca,'FontSize',fs)
                %
                figure_name = sprintf('Figures/%s_TC%d_tf%1.1f_Convg_Order_s%dp%dq%d_m%d_sp_%d.pdf',eqn,TC,tf,S(2*i),P(2*i),Q(2*i),m,spatial_order);
                print(gcf,figure_name,'-dpdf','-r100','-bestfit')
            end
    end
end


