load('~/Desktop/sinews_homog.mat')
%% get exact marginal
V0 = std(Xscell{1}(:,:,1)).^2;
varfun = @(t,i) 2*w_sparse(3)*t+V0(i);
mu0=mean(Xscell{1}(:,:,1));
sol = @(x,i,t) exp(-x.^2./(2*varfun(t,i)))./sqrt(2*pi*varfun(t,i)); 
i=1;
for j=1:length(xs_obs{end})
    exact_sol(:,j) = sol(xs_obs{1}+c(i)*xs_obs{end}(j)+mu0(i),i,xs_obs{end}(j));
end

%% sinews homogenization example
d=20;
tinds = 2:d:length(xs_obs{end});
figure(1);clf
h1=[];h2=[];
for j=1:length(tinds)
    osc_dat = mean(U_osc{1}(:,:,tinds(j)),2);
    osc_dat = osc_dat/norm(osc_dat,1)/mean(diff(xs_osc{1}));
    h1=[h1;fill(xs_osc{1}+2.5*(j-1),osc_dat,[0 0 0],'facealpha',0.50,'linewidth',1)];
    hold on
    % homog_dat = mean(U_obs{1}(:,:,tinds(j)),2);
    homog_dat = exact_sol(:,tinds(j));
    homog_dat = homog_dat/norm(homog_dat,1)/mean(diff(xs_obs{1}));
    h2=[h2;fill(xs_obs{1}+2.5*(j-1),homog_dat,[0 1 1],'facealpha',0.50,'linewidth',1.8)];
end
hold off
xlim([-3.5 12])
ylim([0 0.7])
legend([h1(1);h2(1)],{'Highly Osc. data','Learned Homog. system'},'box','off','interpreter','latex',...
    'location',[0.5 0.75 0.3 0.15],'textcolor','black')
xlabel('$x$','interpreter','latex'); ylabel('Particle density $\rho$','interpreter','latex')
set(gca,'ticklabelinterpreter','latex','fontsize',18)

x = [-1.1];
y = [0.08];
text(x,y,'$t=0.04$','interpreter','latex','fontsize',14,'Rotation',45)
x = [2.5];
y = [0.06];
text(x,y,'$t=0.84$','interpreter','latex','fontsize',14,'Rotation',45)
x = [6.4];
y = [0.05];
text(x,y,'$t=1.64$','interpreter','latex','fontsize',14,'Rotation',45)
set(gca,'position',[0.1315    0.1409    0.77    0.77])
% x = [0.6 0.7];
% y = [0.65 0.65];
% annotation('textarrow',x,y,'String','{\it time}','interpreter','latex')

saveas(gcf,['~/Desktop/homog.png'])
%% gather data -  first run local2D_livescript;
U_osc = U_obs;
xs_osc = xs_obs;

w_sparse = W(W~=0);
sig = w_sparse(3); % diffusion coefficient 
d = 2; % spatial dimension
c = w_sparse(1:2)'; % drift velocity
N = 2^18; % number of particles 
numt = length(xs_obs{end}); % number of timepoints
dt = mean(diff(xs_obs{end})); % timestep
t = 0:dt:(numt-1)*dt; % timegrid
N0 = std(Xscell{1}(:,:,1)).*randn(N,d)+mean(Xscell{1}(:,:,1)); % initial conditions
Xscell = {N0+sqrt(2*dt*sig)*cumsum(randn(N,d,numt),3) - dt*c.*reshape((0:(numt-1)),1,1,numt)};

exps = 1; %choose experiments from Xscell to use in regression (if Xscell contains multiple exps)
NN = N; %choose number of particles to use from each experiment
numx = 200 + 1; %choose histogram grid resolution
numsdv = inf;
custdom = [];
coarsen_data = [[0 1 1];[0 1 1];[0 1 1]];
scoord = 0;
Xsnz = [0 1]; %set extrinsic noise level

get_particle_distrib; %compute histogram


