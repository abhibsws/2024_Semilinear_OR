% Different test problems
EQN = {'SemiLinAdv','VisBurgers','AllenCahn2D'};

% Different initial conditions
IC = [1,1,1];
% Final times
TF = [0.5,1,0.5];
% Spatial order 
sp_or = {5,6,0};
% Number of grid points
M = [2048,1e3,25];

%plot colors and linestyles
C = {'b','r','m'}; Cref = {[0.5,0.5,0.5]};
linS = {'-','--',':'}; Mar = {'o','s','+'}; ms = 10; fs = 20;

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
            S=[8,5,19,5,10,5]; P=[4,4,5,5,5,5]; Q=[3,1,4,1,4,1]; SchNo=[3,1,6,2,5,2];
            filename= sprintf('errdata/%s/ErrConvgData_%s_TC%d_tf%.1f.mat',eqn,eqn,TC,tf);
            load(filename)

            st = [1,4,1,4,1,4]; en = [4,7,4,7,4,7];
            Cof = [1e-1,4e1,4e-1,4e2,6e-1,1e3]; Sl = [2,4,2,5,2,5]; 
            yl = [5e-13,1e-13,1e-13]; yu = [5e-4,1e-3,1e-3];
            % Plot
            for i = 1:3
                figure(i)
                set(gcf,'position',[0 0 600 600])
                set(gca, 'LooseInset', max(get(gca, 'TightInset'), 0)); % remove white space
                set(0,'DefaultLineLineWidth',3);
                loglog(dts(2*i-1,:),U_Err(2*i-1,:),'linestyle',linS{1},'color',C{1},'marker',Mar{1},'MarkerSize',ms)
                hold on
                loglog(dts(2*i,:),U_Err(2*i,:),'linestyle',linS{2},'color',C{1},'marker',Mar{2},'MarkerSize',ms)
                hold on
                loglog(dts(2*i-1,st(2*i-1):en(2*i-1)),Cof(2*i-1)*dts(2*i-1,st(2*i-1):en(2*i-1)).^Sl(2*i-1),'--','color',Cref{1})
                loglog(dts(2*i,st(2*i):en(2*i)),Cof(2*i)*dts(2*i,st(2*i):en(2*i)).^Sl(2*i),'-','color',Cref{1})
                legend('','Location', 'northwest')

                legend({sprintf('$\\mathbf{\\textbf{(%d,%d,%d)}}$', S(2*i-1), P(2*i-1), Q(2*i-1)),...
                sprintf('(%d,%d,%d)', S(2*i), P(2*i), Q(2*i)),...
                sprintf('Slope %d', Sl(2*i-1)),...
                sprintf('Slope %d', Sl(2*i))},...
                'NumColumns', 2, 'Interpreter', 'latex');
                %
                xlim([dts(1,end),dts(1,1)])
                ylim([yl(i),yu(i)])
                xlabel('\Delta t');
                %ylabel('Error');
                ylabel('Error', 'Position', get(gca, 'YLabel').Position + [-1e-5, 1e-7, 0])
                grid minor
                set(gca,'FontSize',fs+5)
                %
                figure_name = sprintf('Figures/%s_TC%d_tf%1.1f_m%d_sp_%d_s%dp%dq%d.pdf',eqn,TC,tf,m,spatial_order,S(2*i-1),P(2*i-1),Q(2*i-1));
                print(gcf,figure_name,'-dpdf','-r100','-vector','-bestfit')
            end
        case 'VisBurgers'
            S=[8,5,10,5]; P=[4,4,5,5]; Q=[3,1,4,1]; SchNo=[3,1,5,2];
            filename= sprintf('errdata/%s/ErrConvgData_%s_TC%d_tf%.1f.mat',eqn,eqn,TC,tf);
            load(filename) 

            st = [1,4,1,4]; en = [4,7,4,7];
            Cof = [2e-1,1e1,2e0,4e2]; Sl = [2,4,2,5]; 
            yl = [2e-12,3e-13]; yu = [3e-3,7e-2]; 
            % Plot
            for i = 1:2
                figure(3+i)
                set(gcf,'position',[0 0 600 600])
                set(gca, 'LooseInset', max(get(gca, 'TightInset'), 0)); % remove white space
                set(0,'DefaultLineLineWidth',3);
                loglog(dts(2*i-1,:),U_Err(2*i-1,:),'linestyle',linS{1},'color',C{1},'marker',Mar{1},'MarkerSize',ms)
                hold on
                loglog(dts(2*i,:),U_Err(2*i,:),'linestyle',linS{2},'color',C{1},'marker',Mar{2},'MarkerSize',ms)
                hold on
                loglog(dts(2*i-1,st(2*i-1):en(2*i-1)),Cof(2*i-1)*dts(2*i-1,st(2*i-1):en(2*i-1)).^Sl(2*i-1),'--','color',Cref{1})
                loglog(dts(2*i,st(2*i):en(2*i)),Cof(2*i)*dts(2*i,st(2*i):en(2*i)).^Sl(2*i),'-','color',Cref{1})
                legend('','Location', 'best')

                legend({sprintf('$\\mathbf{\\textbf{(%d,%d,%d)}}$', S(2*i-1), P(2*i-1), Q(2*i-1)),...
                sprintf('(%d,%d,%d)', S(2*i), P(2*i), Q(2*i)),...
                sprintf('Slope %d', Sl(2*i-1)),...
                sprintf('Slope %d', Sl(2*i))},...
                'NumColumns', 2, 'Interpreter', 'latex');
                %
                xlim([dts(1,end),dts(1,1)])
                ylim([yl(i),yu(i)])
                xlabel('\Delta t');
                ylabel('Error', 'Position', get(gca, 'YLabel').Position + [-1e-6, 5e-8, 0])
                grid minor
                set(gca,'FontSize',fs)
                %
                figure_name = sprintf('Figures/%s_TC%d_tf%1.1f_m%d_sp_%d_s%dp%dq%d.pdf',eqn,TC,tf,m,spatial_order,S(2*i-1),P(2*i-1),Q(2*i-1));
                print(gcf,figure_name,'-dpdf','-r100','-vector','-bestfit')
            end
        case 'AllenCahn2D'
            S=[8,5,10,5]; P=[4,4,5,5]; Q=[3,1,4,1]; SchNo=[3,1,5,2];
            filename= sprintf('errdata/%s/ErrConvgData_%s_TC%d_tf%.1f.mat',eqn,eqn,TC,tf);
            load(filename)

            st = [1,4,1,4]; en = [4,7,4,7];
            Cof = [1e0,1e2,1e1,1e3]; Sl = [2,4,2,5]; 
            yl = [1e-11,1e-13]; yu = [2e-2,1e-1];
            yl_pos = [-1.9e-7,5e-8];
            % Plot
            for i = 1:2
                figure(5+i)
                set(gcf,'position',[0 0 600 600])
                set(gca, 'LooseInset', max(get(gca, 'TightInset'), 0)); % remove white space
                set(0,'DefaultLineLineWidth',3);
                loglog(dts(2*i-1,:),U_Err(2*i-1,:),'linestyle',linS{1},'color',C{1},'marker',Mar{1},'MarkerSize',ms)
                hold on
                loglog(dts(2*i,:),U_Err(2*i,:),'linestyle',linS{2},'color',C{1},'marker',Mar{2},'MarkerSize',ms)
                hold on
                loglog(dts(2*i-1,st(2*i-1):en(2*i-1)),Cof(2*i-1)*dts(2*i-1,st(2*i-1):en(2*i-1)).^Sl(2*i-1),'--','color',Cref{1})
                loglog(dts(2*i,st(2*i):en(2*i)),Cof(2*i)*dts(2*i,st(2*i):en(2*i)).^Sl(2*i),'-','color',Cref{1})
                legend('','Location', 'best')

                legend({sprintf('$\\mathbf{\\textbf{(%d,%d,%d)}}$', S(2*i-1), P(2*i-1), Q(2*i-1)),...
                sprintf('(%d,%d,%d)', S(2*i), P(2*i), Q(2*i)),...
                sprintf('Slope %d', Sl(2*i-1)),...
                sprintf('Slope %d', Sl(2*i))},...
                'NumColumns', 2, 'Interpreter', 'latex');
                %
                xlim([dts(1,end),dts(1,1)])
                ylim([yl(i),yu(i)])
                xlabel('\Delta t');
                ylabel('Error', 'Position', get(gca, 'YLabel').Position + [-1e-6, yl_pos(i), 0])
                grid minor
                set(gca,'FontSize',fs)

                figure_name = sprintf('Figures/%s_TC%d_tf%1.1f_m%d_sp_%d_s%dp%dq%d.pdf',eqn,TC,tf,m,spatial_order,S(2*i-1),P(2*i-1),Q(2*i-1));
                print(gcf,figure_name,'-dpdf','-r100','-vector','-bestfit')
            end
    end
end


