eq_ind = 3; % change here for different test cases
% Different test problems
EQN = {'SemiLinAdv','VisBurgers','AllenCahn2D'};

% Initial conditions
IC = [1,1,1];
% Final times
TF = [0.5,1,0.5];
% Spatial orders
sp_or = {5,6,'Cheb'};
% Number of grid points
M = [5e3,1e3,30];

eqn=EQN{eq_ind};
TC=IC(eq_ind); tf=TF(eq_ind); m = M(eq_ind); spatial_order = sp_or{eq_ind};

switch eqn 
    case 'SemiLinAdv'
        % DIRK-(s,p,q) scheme
        S=[5,7,5,10]; 
        P=[4,4,5,5]; 
        Q=[1,4,1,4]; 
        SchNo=[2,9,3,10];
        NT = 2.^(4:1:10); %number of time steps dt=1/NT
    case 'VisBurgers'
        % DIRK-(s,p,q) scheme
        S=[5,7,7,5,12,14]; 
        P=[4,4,4,5,5,5]; 
        Q=[1,4,4,1,5,4]; 
        SchNo=[2,4,9,3,6,11];
        NT = 2.^(3:1:10);
    case 'AllenCahn2D'
        % DIRK-(s,p,q) scheme
        S=[5,7,5,10]; 
        P=[4,4,5,5]; 
        Q=[1,4,1,4]; 
        SchNo=[2,9,3,10];
        NT = 2.^(3:1:10);
end
% Compute
dts = zeros(length(S),length(NT)); 
U_Err = zeros(length(S),length(NT)); dU_Err = zeros(length(S),length(NT));

% Create a folder for each test case in 'errdata' folder
save_data = 1;
if save_data
    foldername_err = sprintf('errdata/%s/',eqn);
    if exist(foldername_err,'dir')==0
        mkdir(foldername_err);
    end
end

tic
for i = 1:length(S)
    s=S(i); p=P(i); q=Q(i); scheme_no=SchNo(i); 
    disp(s)
    parfor j = 1:length(NT)
        dt = (1/NT(j)); nt = ceil(tf/dt);dt = tf/nt;
        display(nt)
        switch eqn 
            case 'SemiLinAdv'
                [uerr,duerr] = DIRK_SemiLinAdvEqn(TC,tf,m,spatial_order,s,p,q,scheme_no,dt,nt);
                dts(i,j) = dt; U_Err(i,j) = uerr; dU_Err(i,j) = duerr;                
            case 'VisBurgers'
                [uerr,duerr] = DIRKspqVisBurgerEqn(TC,tf,m,spatial_order,s,p,q,scheme_no,dt,nt);
                dts(i,j) = dt; U_Err(i,j) = uerr; dU_Err(i,j) = duerr; 
            case 'AllenCahn2D'
              uerr = DIRK_Allen_Cahn_2D_Eqn(TC,tf,m,spatial_order,s,p,q,scheme_no,dt,nt);
                dts(i,j) = dt; U_Err(i,j) = uerr;
        end
    end
end
toc
%-------------------------------------------------------------------------%
% Saving data 
if save_data
    filename= sprintf('%sErrConvgData_%s_TC%d_tf%.1f.mat',...
        foldername_err,eqn,TC,tf);
    switch eqn 
        case 'SemiLinAdv'
            save(filename, 'dts', 'U_Err', 'dU_Err');
        case 'VisBurgers'
            save(filename, 'dts', 'U_Err', 'dU_Err');
        case 'AllenCahn2D'
            save(filename, 'dts', 'U_Err');
    end
end
%-------------------------------------------------------------------------%
