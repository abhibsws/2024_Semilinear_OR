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
for eq_ind = 3:3
    eqn=EQN{eq_ind};TC=IC(eq_ind); tf=TF(eq_ind); 
    m = M(eq_ind); spatial_order = sp_or{eq_ind};
    % for legends
    switch eqn 
        case 'SemiLinAdv'
            S=[8,5,19,5,10,5]; P=[4,4,5,5,5,5]; Q=[3,1,4,1,4,1]; SchNo=[3,1,6,2,5,2];
            filename= sprintf('errdata/%s/ErrConvgData_%s_TC%d_tf%.1f.mat',eqn,eqn,TC,tf);
            load(filename)

            st = [1,4,1,4,1,4]; en = [4,7,4,7,4,7];
            Cof = [1e-1,4e1,4e-1,4e2,4e-1,4e2]; Sl = [2,4,2,5,2,5]; 
            yl = [1e-13,1e-13,1e-13]; yu = [1e-3,1e-3,1e-3];
            % Coordinates to modify the positions of the xlabel and ylabel
            mod_XL = [
                      [0, 8e-14, 0];
                      [0, 2e-14, 0];
                      [0, 2e-14, 0]
                      ];
            mod_YL = [
                      [-1e-5, 1e-7, 0];
                      [-1e-5, 1e-7, 0];
                      [-1e-5, 1e-7, 0]
                      ];
            % Plot
            for i = 1:3
                %---------for minimum white space
                % Ensure figure remains valid
                fig = figure(i);
                % Set figure size
                set(fig, 'Units', 'pixels', 'Position', [50, 50, 600, 600]);
                % Adjust axis to remove extra margins
                ax = gca;
                set(ax, 'LooseInset', [0, 0, 0, 0]);
                set(0,'DefaultLineLineWidth',3);
                %---------for minimum white space
                loglog(dts(2*i-1,:),U_Err(2*i-1,:),'linestyle',linS{1},'color',C{1},'marker',Mar{1},'MarkerSize',ms+2)
                hold on
                loglog(dts(2*i,:),U_Err(2*i,:),'linestyle',linS{2},'color',C{1},'marker',Mar{2},'MarkerSize',ms+2)
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
                xlabel('h', 'Position', get(gca, 'XLabel').Position + mod_XL(i,:))
                ylabel('Error', 'Position', get(gca, 'YLabel').Position + mod_YL(i,:))
                grid minor
                set(gca,'FontSize',fs+5)
                %
                figure_name = sprintf('Figures/%s_TC%d_tf%1.1f_m%d_sp_%d_s%dp%dq%d.pdf',eqn,TC,tf,m,spatial_order,S(2*i-1),P(2*i-1),Q(2*i-1));
                % print(gcf,figure_name,'-dpdf','-r100','-vector','-bestfit')
                exportgraphics(gcf, figure_name, 'ContentType', 'vector')
            end
        case 'VisBurgers'
            S=[8,5,10,5]; P=[4,4,5,5]; Q=[3,1,4,1]; SchNo=[3,1,5,2];
            filename= sprintf('errdata/%s/ErrConvgData_%s_TC%d_tf%.1f.mat',eqn,eqn,TC,tf);
            load(filename) 

            st = [1,4,1,4]; en = [4,7,4,7];
            Cof = [2e-1,1e1,2e0,4e1]; Sl = [2,4,2,4.5]; 
            yl = [2e-12,3e-13]; yu = [3e-3,1e-1]; 
            % Coordinates to modify the positions of the xlabel and ylabel
            mod_XL = [
                      [-5e-3, 4e-13, 0];
                      [-5e-3, 5e-14, 0]
                      ];
            mod_YL = [
                      [1e-4, 5e-8, 0];
                      [1e-4, -0.3e-7, 0]
                      ];
            % Plot
            for i = 1:2
                %---------for minimum white space
                % Ensure figure remains valid
                fig = figure(3+i);
                % Set figure size
                set(fig, 'Units', 'pixels', 'Position', [50, 50, 500, 500]);
                % Adjust axis to remove extra margins
                ax = gca;
                set(ax, 'LooseInset', [0, 0, 0, 0]);
                set(0,'DefaultLineLineWidth',3);
                %---------for minimum white space
                loglog(dts(2*i-1,:),U_Err(2*i-1,:),'linestyle',linS{1},'color',C{1},'marker',Mar{1},'MarkerSize',ms)
                hold on
                loglog(dts(2*i,:),U_Err(2*i,:),'linestyle',linS{2},'color',C{1},'marker',Mar{2},'MarkerSize',ms)
                hold on
                loglog(dts(2*i-1,st(2*i-1):en(2*i-1)),Cof(2*i-1)*dts(2*i-1,st(2*i-1):en(2*i-1)).^Sl(2*i-1),'--','color',Cref{1})
                loglog(dts(2*i,st(2*i):en(2*i)),Cof(2*i)*dts(2*i,st(2*i):en(2*i)).^Sl(2*i),'-','color',Cref{1})
                legend('','Location', 'best')
                %
                if i == 1
                    legend({sprintf('$\\mathbf{\\textbf{(%d,%d,%d)}}$', S(2*i-1), P(2*i-1), Q(2*i-1)),...
                    sprintf('(%d,%d,%d)', S(2*i), P(2*i), Q(2*i)),...
                    sprintf('Slope %d', Sl(2*i-1)),...
                    sprintf('Slope %d', Sl(2*i))},...
                    'NumColumns', 2, 'Interpreter', 'latex');
                elseif i == 2
                    legend({sprintf('$\\mathbf{\\textbf{(%d,%d,%d)}}$', S(2*i-1), P(2*i-1), Q(2*i-1)),...
                    sprintf('(%d,%d,%d)', S(2*i), P(2*i), Q(2*i)),...
                    sprintf('Slope %d', Sl(2*i-1)),...
                    sprintf('Slope %1.1f', Sl(2*i))},...
                    'NumColumns', 2, 'Interpreter', 'latex');
                end
                %
                xlim([dts(1,end),dts(1,1)])
                ylim([yl(i),yu(i)])
                xlabel('h', 'Position', get(gca, 'XLabel').Position + mod_XL(i,:))
                ylabel('Error', 'Position', get(gca, 'YLabel').Position + mod_YL(i,:))
                grid minor
                set(gca,'FontSize',fs)
                %
                figure_name = sprintf('Figures/%s_TC%d_tf%1.1f_m%d_sp_%d_s%dp%dq%d.pdf',eqn,TC,tf,m,spatial_order,S(2*i-1),P(2*i-1),Q(2*i-1));
                % print(gcf,figure_name,'-dpdf','-r100','-vector','-bestfit')
                exportgraphics(gcf, figure_name, 'ContentType', 'vector')
                %--------- compute slopes
                xdata = dts(2*i-1, :);  % Grid sizes (h)
                ydata = U_Err(2*i-1, :); % Errors
                
                logx = log(xdata);
                logy = log(ydata);
            
                slopes = diff(logy) ./ diff(logx); % Vector of local slopes
                VisBurgers_computed_slopes{i} = slopes; % Store it in the cell array
            end
        case 'AllenCahn2D'
            S=[8,5,10,5]; P=[4,4,5,5]; Q=[3,1,4,1]; SchNo=[3,1,5,2];
            filename= sprintf('errdata/%s/ErrConvgData_%s_TC%d_tf%.1f.mat',eqn,eqn,TC,tf);
            load(filename)

            st = [4,4,1,4]; en = [6,8,4,8];
            Cof = [8e-1,5e1,1e1,1e3]; Sl = [2,4,2,5]; 
            yl = [1e-11,1e-13]; yu = [1e-2,1e-1];
            % Coordinates to modify the positions of the xlabel and ylabel
            mod_XL = [
                      [-5e-3, 2e-12, 0];
                      [-5e-3, 2e-14, 0]
                      ];
            mod_YL = [
                      [1e-4, -2e-7, 0];
                      [1e-4, 0, 0]
                      ];
            % Plot
            for i = 1:2
                %---------for minimum white space
                % Ensure figure remains valid
                fig = figure(5+i);
                % Set figure size
                set(fig, 'Units', 'pixels', 'Position', [50, 50, 500, 500]);
                % Adjust axis to remove extra margins
                ax = gca;
                set(ax, 'LooseInset', [0, 0, 0, 0]);
                set(0,'DefaultLineLineWidth',3);
                %---------for minimum white space
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
                xlabel('h', 'Position', get(gca, 'XLabel').Position + mod_XL(i,:))
                ylabel('Error', 'Position', get(gca, 'YLabel').Position + mod_YL(i,:))
                grid minor
                set(gca,'FontSize',fs)

                figure_name = sprintf('Figures/%s_TC%d_tf%1.1f_m%d_sp_%d_s%dp%dq%d.pdf',eqn,TC,tf,m,spatial_order,S(2*i-1),P(2*i-1),Q(2*i-1));
                % print(gcf,figure_name,'-dpdf','-r100','-vector','-bestfit')
                exportgraphics(gcf, figure_name, 'ContentType', 'vector')
            end
    end
end


