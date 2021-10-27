%% particle data - Xscell, t

exps = 1:8; 
subsampN = 1;

numx = 256 + 1;
numsdv = 3;
custdom = [];
coarsen_data = [[0 1 1];[0 1 1];[0 1 1]];
Xsnz = [0 1];
scoord = 0;

NN = size(Xscell{1},1);
dim = size(Xscell{1},2)+1;

%% Set Library

% local autonomous operators
max_dx = 0;
max_dt = 0;
polys = [];
trigs = [];
use_all_dt = 0;
use_cross_dx = 0;
custom_add = [];
custom_remove = {@(mat,lhs) find(all([mat(:,2)==0 mat(:,3)==0 ~ismember(mat,lhs,'rows')],2))};
toggle_comb = 1;

% local non-autonomous operators
driftpolys=[0 2:8]; drifttrigs=[]; 
diffpolys=[]; difftrigs=[0:7]; crossdrift=0;

% non-local operators
convargs = {{'dimx',dim-1,'utagin',1,'utagout',1,'psitags',1,'Mon',[0:7],'Sing',[],'Exp',[],'Singeps',[],'svdtol',[1e-8]}};

%% Set Discretization

phi_class = {1,1};
sm_x = 4;
sm_t = 4;

% %%% uncomment to set test fcn params manually
% p_x = 6; mxs = 12;
% tau_x = []; k_x = []; tauhat_x = []; 
% p_t = 4; mts = 4;
% tau_t = []; k_t = []; tauhat_t = [];

%%% set test fcn params using cornerpoint
tauhat_x = 1; tau_x = 10^-4;
p_x = []; mxs = 1;
tauhat_t = 1; tau_t = 10^-4; 
p_t = []; mts = 1;

%%% rescale coordinates and/or use approx. variance for improved conditioning
scales = 0;
covtol = 0;

%%% trim rows with low particle density
trim_tags = [1 zeros(1,dim-1) 0];
trim_fcn = {@(col) max(col/max(abs(col)),eps)};
inds_keep_fcn = @(col) log10(col)>-1;

%% MSTLSQP

lambda = 10.^(linspace(-4, -0, 100));
gamma_tol = Inf;
alpha = 0.5;
maxits = 10;
sparsity_scale = 0;
excl_tags={'u^{1}_{lap}','u^{1}_{xx}'};
excl_tols=0;
tol = 10^-8; 
maxQPits = 100; 
dispQP = 'off';
meth = 'STLSQP';

%% run alg

wsindy_IPS_script;

%% Display 

print_loc = 1;
toggle_plot_basis_fcn = 0;
toggle_plot_sol = 0;
toggle_plot_loss = 0;
toggle_plot_drift = 0;
toggle_plot_IPforce = 0;
toggle_plot_fft = {[1 2 3],1,1};

get_results;
wsindy_IPS_display;