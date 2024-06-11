%% spatial domain
Nx = 80;
x = linspace(xs_obs{1}(1),xs_obs{1}(end),Nx);
dx = mean(diff(x));
y = xs_obs{2}(1):dx:xs_obs{2}(end);
[xx,yy] = meshgrid(x,y);%%% x and y are flipped
Ny = length(y); 

%% initial conditions
[xx_0,yy_0] = meshgrid(xs_obs{1},xs_obs{2});
rho0 = interp2(xx_0,yy_0,U_obs{1}(:,:,1)',xx,yy);
rho0 = smoothdata(rho0,1,"gaussian",Ny/20);
rho0 = smoothdata(rho0,2,"gaussian",Nx/20);
% l = 10;p=2;
% rho0(1:l,:) = rho0(1:l,:).*((1-(l-1:-1:0)/(l-1)).^p)';
% rho0(end-l+1:end,:) = rho0(end-l+1:end,:).*((1-(l-1:-1:0)/(l-1)).^p)';
figure(1);
surf(xx,yy,rho0)
view(0, 90)
xlabel x, ylabel y
axis equal

%% local forces
% nu_mat = 0.02+0*xx; %fzero(@ (ep) 4*dx^2-ep*(-2*log(eps)-log(2*pi*ep)),1);
nu_mat = 10*sigf(xx,yy,0);
[D1x,D2x] = getDmats(x);
[D1y,D2y] = getDmats(y);

Ex = speye(Nx); Ex(1,1) = 0; Ex(end,end) = 0;
Ey = speye(Ny); Ey(1,1) = 0; Ey(end,end) = 0;
Dx = kron(D1x,Ey); Dy = kron(Ex,D1y);
Dlap = kron(Ex,D2y)+kron(D2x,Ey);
Dlap_lin = Dlap*spdiags(nu_mat(:),0,Nx*Ny);

%% growth term
fg = @(u) u*0;
w = W;
pows = [];
coeffs = [];
lin_growth = 0;
for j=1:length(tags_pde_G)
    if isequal(tags_pde_G{j}(end-2:end),'_{}')
        if w(j)~=0
            pow = str2num(tags_pde_G{j}(end-4));
            if pow~=1
                pows = [pows pow];
                coeffs = [coeffs w(j)];
                fg = @(u) fg(u)+w(j)*u.^pows(end);
            else
                lin_growth = w(j);
            end
        end
    end
end
growthmat = ones(Ny,Nx);
growthmat(:,1) = 0; growthmat(:,end) = 0;
growthmat(1,:) = 0; growthmat(end,:) = 0;
growthmat = spdiags(growthmat(:),0,Nx*Ny,Nx*Ny);

figure(2)
uplot=linspace(min(rho0(:)),2*max(rho0(:)),100);
plot(uplot,fg(uplot)+lin_growth*uplot)

%% nonlocal forces
fx = @(x,y)f_learnedcell{1}(x,y,0);
fy = @(x,y)f_learnedcell{2}(x,y,0);

zz_x = fliplr([x(end)-x x(1)-x(2:end)]);
zz_y = fliplr([y(end)-y y(1)-y(2:end)]);
[xxx,yyy] = meshgrid(zz_x,zz_y);
rad = hypot(xxx,yyy);
forceX = fx(xxx,yyy);
forceY = fy(xxx,yyy);

figure(4)
subplot(1,3,1)
plot(forceX(floor(end/2),:)'.*[1 0])
subplot(1,3,2)
plot(forceY(:,floor(end/2))*[1 0])
subplot(1,3,3)
imagesc(zz_x,zz_y,hypot(forceX,forceY))
axis equal
colorbar

%% set bcs

% bc='Dirichlet';
bc='Neumann';

%% sim

gap = 20;
dt = mean(diff(xs_obs{end}))/gap;
M = length(xs_obs{end});

t = 0:dt*gap:(M-1)*dt*gap;
rhos = zeros(Ny,Nx,M);
rho = rho0;
rhos(:,:,1) = rho;

A0 = kron(Ex,Ey) - (dt*lin_growth)*growthmat - dt*Dlap_lin;
BC_mask = sparse(Ny,Nx);
BC_mask(:,1) = 1;
BC_mask(:,end) = 1;
BC_mask(1,:) = 1;
BC_mask(end,:) = 1;
NLrad = 1;
nl_mask = (1-((NLrad:-1:0)/NLrad).^2).^2;

for tt = 1:(M-1)*gap
    tic
    vx = conv2(rho,forceX,'same')*dx^2;
    vy = conv2(rho,forceY,'same')*dx^2;

    vx(:,1:NLrad+1) = vx(:,1:NLrad+1).*nl_mask;
    vx(:,end-NLrad:end) = vx(:,end-NLrad:end).*fliplr(nl_mask);
    vy(1:NLrad+1,:) = vy(1:NLrad+1,:).*nl_mask';
    vy(end-NLrad:end,:) = vy(end-NLrad:end,:).*fliplr(nl_mask)';

    if isequal(bc,'Neumann')
        BCmat_x_L = nu_mat(:,1:2)/dx.*[-1 1] + vx(:,1:2)/2;
        BCmat_x_R = nu_mat(:,end-1:end)/dx.*[-1 1] + vx(:,end-1:end)/2;
        BC_temp = sparse(Nx,Nx);
        BC_temp(1,1:2)=1;
        BC_temp(end,end-1:end)=1;
        BCmat_x = kron(BC_temp,speye(Ny));
        BCmat_x(find(BCmat_x)) = [BCmat_x_L';BCmat_x_R']';

        BCmat_y_D = nu_mat(1:2,:)/dx.*[-1;1] + vy(1:2,:)/2;
        BCmat_y_U = nu_mat(end-1:end,:)/dx.*[-1;1] + vy(end-1:end,:)/2;
        BC_temp = sparse(Ny,Ny);
        BC_temp(1,1:2)=1;
        BC_temp(end,end-1:end)=1;
        BCmat_y = kron(speye(Nx),BC_temp);
        BCmat_y(find(BCmat_y)) = [BCmat_y_D;BCmat_y_U];

        BCmat = BCmat_x+BCmat_y;
        BC_rhs = zeros(size(BCmat,1),1);

        % BC_rhs = rho.*BC_mask;

    elseif isequal(bc,'Dirichlet')
        BC_temp = sparse(Nx,Nx);
        BC_temp(1,1)=1;
        BC_temp(end,end)=1;
        BCmat_x = kron(BC_temp,speye(Ny));

        BC_temp = sparse(Ny,Ny);
        BC_temp(1,1)=1;
        BC_temp(end,end)=1;
        BCmat_y = kron(speye(Nx),BC_temp);

        BCmat = BCmat_x+BCmat_y;
        BCmat(1,1) = 1;
        BCmat(Ny,Ny) = 1;
        BCmat(end,end) = 1;
        BCmat((Nx-1)*Ny+1,(Nx-1)*Ny+1) = 1;

        % BC_rhs = interp2(xx_0,yy_0,U_obs{1}(:,:,ceil(tt/gap))',xx,yy)*(ceil(tt/gap)-tt/gap) ...
        %     + interp2(xx_0,yy_0,U_obs{1}(:,:,ceil(tt/gap)+1)',xx,yy)*(tt/gap-floor(tt/gap));
        BC_rhs = interp2(xx_0,yy_0,U_obs{1}(:,:,ceil(tt/gap))',xx,yy);
        BC_rhs = BC_rhs.*BC_mask;
        BC_rhs = BC_rhs(:);
    end

    growth_term = fg(rho(:));
    flux = Dx*spdiags(vx(:), 0, Nx*Ny, Nx*Ny)+Dy*spdiags(vy(:), 0, Nx*Ny, Nx*Ny);
    fg_t = fg(rho);
    fg_t(:,1) = 0;
    fg_t(:,end) = 0;
    fg_t(1,:) = 0;
    fg_t(end,:) = 0;
    rho(:,1) = 0;
    rho(:,end) = 0;
    rho(1,:) = 0;
    rho(end,:) = 0;
    RHS = rho(:) + dt*fg_t(:) + BC_rhs;
    A = A0 - dt*flux + BCmat;

    % rho = A\RHS;
    [rho,fl] = gmres(A,RHS,[],[],[],[],[],rho(:));
    rho = reshape(rho,Ny,Nx);
    disp([tt/M/gap toc max(rho(:))])
    if mod(tt,gap)==0 
        rhos(:,:,tt/gap+1) = rho;
        plot(y,mean(rho,2));
        hold on
        plot(xs_obs{2},mean(U_obs{1}(:,:,tt/gap)',2))
        hold off
        ylim([0 2*max(rho0(:))])
        drawnow
    end
end

%% put in WSINDy form
U_exact = {rhos};
xs = {y,x,t};
lhs = [1 0 0 1];

%% view 2D
figure(8)
for j=1:length(t)
    subplot(2,1,1)
    surf(x,y,U_exact{1}(:,:,j),'edgecolor','none')
    view([0 90])
    % xlim([50 350])
    % ylim([100 950])
    zlim([0 2])
    subplot(2,1,2)
    surf(xs_obs{1},xs_obs{2},U_obs{1}(:,:,j)','edgecolor','none')
    view([0 90])
    % xlim([50 350])
    % ylim([100 950])
    zlim([0 2])
    % axis equal
    % colorbar
    drawnow
end

%% view 1D

for j=1:length(t)
    plot(y,mean(U_exact{1}(:,:,j),2));
    hold on
    plot(xs_obs{2},mean(U_obs{1}(:,:,j)',2))
    hold off
    ylim([0 2])
    % pause(0.1)
    drawnow
end

function [D1,D2] = getDmats(x)
    N = length(x);
    dx = mean(diff(x));
    e = ones(N,1);     
    D1 = spdiags([-e zeros(N,1) e], -1:1, N, N);
    D2 = spdiags([e -2*e e], -1:1, N, N);
    D1(1,1:2) = 0; D2(1,1:2) = 0;
    D1(N,N-1:N) = 0; D2(N,N-1:N) = 0;
    D1 = D1/2/dx;
    D2 = D2/dx^2;
end
% vx = -conv2(rho,forceX,'same')*dx^2;
% vy = -conv2(rho,forceY,'same')*dx^2;
% dt = 10*sqrt(2)*dx/max(max(hypot(vx,vy)));
% cfl = sqrt(2)*dx/dt;

