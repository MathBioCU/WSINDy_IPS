d = 2;
nu = 1/8/pi;
N = 100;
dt = 0.05;
M = 400;
subdt = 20;
tol_ode = 10^-8;
verbose = 1;
bc = -1;

rng('shuffle')

mus=[1 1];
sigs = zeros(size(mus,2)*size(mus,1),size(mus,2));
for i=1:size(mus,1)
    sigtemp = 2*eye(d);
    sigs([size(mus,2)*(i-1)+1:size(mus,2)*i],:) = sigtemp;
end
x0_dist = 2;

delt = 0.01; reg = 0; p = 0; scal = 1;
[K_true,f_true] = QANR_reg(delt,reg,p,scal,d);
f_true = @(x,t)f_true(x); %%% assumed radial component

% positional force
g_true = @(x,t) [-0.05 -0.025]; %%% assumed vector component
g_true = [];

% vortical force
thet = 0;
A = @(x,y) (exp(1i*thet)*(x+1i*y)).^4/2;
Ax = @(x,y) real(A(x,y));
Ay = @(x,y) -imag(A(x,y));
Ax = [];
Ay = [];

RHS = @(X,t) rhs(f_true,g_true,Ax,Ay,X,d,t);
sig = @(X,t) X*0+sqrt(2*nu);
% sig = @(X,t) sqrt(2*nu)*(2-X);
t = 0:dt:M*dt;

Xtemp = gen_X0(mus,sigs,N,x0_dist);

if nu>=0
    [~,Xtemp,B] = simnu(RHS,t,Xtemp(:),subdt,sig,verbose);
else
    options = odeset('RelTol',tol_ode,'AbsTol',tol_ode*ones(1,length(Xtemp(:))));
    [~,Xtemp]=ode45(@(t,x)RHS(x,t),t,Xtemp(:),options);
end

Xscell = zeros(N,d,length(t));
for m=1:length(t)
X = reshape(Xtemp(m,:)',N,d);
Xscell(:,:,m) = X;
end
Xscell = {Xscell}; 

% save([savedir,'log2D_N',num2str(N),'_nu',num2str(log10(nu)),'_run',FID,'.mat'])


% Functions
function [t,X,B] = simnu(RHS,t,Xtemp,m,sig,verbose)
    M = length(t);
    dt = t(2)-t(1);
    Ntot = length(Xtemp);
    B = zeros(Ntot,M);
    X = zeros(M,Ntot);
    X(1,:) = Xtemp;
    for mm=2:M
        tic,
        if m>1
            t_temp = t(mm-1)+dt/m:dt/m:t(mm);
        else
            t_temp = t(mm-1)+dt;
        end
        for j=1:length(t_temp)
            if ~isempty(sig)
                Btemp = sqrt(dt/m)*sig(Xtemp,t_temp(j)).*randn(Ntot,1);
                B(:,mm) = B(:,mm)+Btemp;
            else
                Btemp = 0*Xtemp;
            end
            Xtemp = Xtemp+dt/m*RHS(Xtemp,t_temp(j))+Btemp;
        end
        X(mm,:) = Xtemp;
        if verbose
            disp([mm toc]);
        end
    end
end    

function [t,X,B] = simnu_hp(RHS,t,Xtemp,m,nu,bc)
    M = length(t);
    dt = t(2)-t(1);
    Ntot = length(Xtemp);
    B = zeros(Ntot,M);
    X = zeros(M,Ntot);
    X(1,:) = Xtemp;
  for mm=2:M
        t_temp = t(mm-1)+dt/m:dt/m:t(mm);
        for j=1:length(t_temp)
            Xtemp = Xtemp+dt/m*RHS(Xtemp,t_temp(j));
            if nu>0
                Btemp = sqrt(2*nu*(dt/m))*randn(Ntot,1);
                Xtemp = Xtemp + Btemp;
                B(:,mm) = B(:,mm)+Btemp;
            end
            if bc == 0
                Xtemp = reshape(Xtemp,[],2);
                Xtemp(:,1) = max(Xtemp(:,1),0);
                Xtemp = Xtemp(:);
            elseif bc == 1
                Xtemp = reshape(Xtemp,[],2);
                Xtemp(:,1) = max(Xtemp(:,1),-Xtemp(:,1));
                Xtemp = Xtemp(:);
            end
        end
        X(mm,:) = Xtemp;
    end
end    

function dX = rhs(f_true,g_true,Ax,Ay,Xtemp,d,t)
    Xtemp = reshape(Xtemp,[],d);
    if ~isempty(f_true)
        update1 = nonlocal_vel(Xtemp,@(x)f_true(x,t));
    else
        update1 = 0*Xtemp;
    end
    if ~isempty(g_true)
        update2 = local_vel(Xtemp,@(x)g_true(x,t));
    else
        update2 = 0*Xtemp;
    end
    if ~isempty(Ax)
        update3 = swirl_vel(Xtemp,Ax,Ay);
    else
        update3 = 0*Xtemp;
    end
    dX = reshape(update1 + update2 + update3,[],1);
end

function dX = nonlocal_vel(agents,f)
        [N,d] = size(agents);
        forces = f(dist(agents'));
        dX = reshape(agents,N,1,d) - reshape(agents,1,N,d);
        dX = dX./max(squeeze(vecnorm(dX,2,3)),10^-16);
       dX = reshape(mean(dX.*forces,2),N,d);
end

function dX = swirl_vel(agents,fx,fy)
        [N,d] = size(agents);
        dX = reshape(agents,N,1,d) - reshape(agents,1,N,d);
        dX = [mean(fx(dX(:,:,1),dX(:,:,2)),2) mean(fy(dX(:,:,1),dX(:,:,2)),2)];
end

function dX = local_vel(agents,V)
        dX = agents*0;
        [N,d] = size(agents);
        for n=1:N
            dX(n,:) = V(agents(n,:));
        end
end

function X0 = gen_X0(mus,sigs,N,x0_dist);
    [m,dim] = size(mus);
    if length(sigs(:))==1;
        sigs = repmat(sigs*eye(dim),m,1);
    end
    X0 = zeros(N,dim);
    inc = floor(N/m);
    for mm=1:m-1             
        if x0_dist == 0
            X0( (mm-1)*inc+(1:inc),:) = mvnrnd(mus(mm,:),sigs(1+(mm-1)*dim:mm*dim,:),inc);
        elseif x0_dist == 1
            X0( (mm-1)*inc+(1:inc),:) = (2*rand(inc,dim)-1)*sigs(1+(mm-1)*dim:mm*dim,:)+mus(mm,:);
        elseif x0_dist == 2
            X0( (mm-1)*inc+(1:inc),:) = samplellipse(inc,sigs(1+(mm-1)*dim,1),sigs(mm*dim,2),mus(mm,1),mus(mm,2));
        end
    end
    inds = (m-1)*inc+1:N;
    if x0_dist == 0
        X0(inds,:) = mvnrnd(mus(m,:),sigs((m-1)*dim+1:end,:),length(inds));
    elseif x0_dist == 1
        X0(inds,:) = (2*rand(length(inds),dim)-1)*sigs(1+(m-1)*dim:end,:)+mus(m,:);
    elseif x0_dist == 2
        X0(inds,:) = samplellipse(length(inds),sigs(1+(m-1)*dim,1),sigs(end,2),mus(m,1),mus(m,2));
    end
end

function X = samplellipse(N,rx,ry,ax,ay)
    init_r = rand(N,1); 
    init_theta = 2*pi*rand(N,1); 
    X = [rx*sqrt(init_r).*cos(init_theta)+ax ry*sqrt(init_r).*sin(init_theta)+ay];
end

function [K,f] = QANR_reg(delta,reg,p,scal,d)

    if ~exist('scal','var')
        scal = 1;
    end
    
    if d==1    
        K = @(x) -scal*abs(x);
        f = @(x) -scal*sign(x);

        if reg > 0
            if reg==1
                K_reg = @(x) -scal*(1/2*delta+1/2/delta*x.^2);
                f_reg = @(x) -scal*(x/delta);
            elseif reg == 2
                K_reg = @(x) -scal*(3/8*delta+3/4/delta*x.^2-1/(8*delta^3)*x.^4);
                f_reg = @(x) -scal*(3/2/delta*x-1/(2*delta^3)*x.^3);
            end

            if mod(p,2) == 0
                K = @(x) K_reg(x).*(abs(x)<delta)+K(x).*(abs(x)>=delta)+x.^p/p;
                f = @(x) f_reg(x).*(abs(x)<delta)+f(x).*(abs(x)>=delta)+x.^(p-1);
            else
                K = @(x) K_reg(x).*(abs(x)<delta)+K(x).*(abs(x)>=delta)+abs(x).^p/p;
                f = @(x) f_reg(x).*(abs(x)<delta)+f(x).*(abs(x)>=delta)+sign(x).*abs(x).^(p-1);
            end

        else
            if mod(p,2) == 0
                K = @(x) K(x) + x.^p/p;
                f = @(x) f(x) + x.^(p-1);
            else
                K = @(x) K(x) + abs(x).^p/p;
                f = @(x) f(x) + sign(x).*abs(x).^(p-1);
            end
        end
    elseif d==2
            
        K = @(x) -1/2/pi*scal*log(max(abs(x),eps));
        f = @(x) -1/2/pi*scal./max(abs(x).^2,eps).*x;
        fp = @(x) 1/2/pi./max(abs(x),eps).^2; 

if and(reg > -inf,delta>0)
        
            if reg<0
                reg =  -reg;
                c_2 = fp(delta)/delta^(reg-2)/reg/(reg-1);
                c_1 = f(delta) - c_2*reg*delta^(reg-1);
                c_0 = K(delta) - c_1*delta - c_2*delta^reg;
                K_reg = @(x) c_0+c_1*x+c_2*max(abs(x),eps).^reg;
                f_reg = @(x) c_1 + c_2*reg*max(abs(x),eps).^(reg-1).*sign(x);
            elseif reg == 0
                c_0 = K(delta)-delta*f(delta);
                c_1 = f(delta);
                K_reg = @(x) c_0+c_1*abs(x);
                f_reg = @(x) c_1*sign(x);    
            elseif reg==1
                c_0 = -(-1+2*log(delta))/(4*pi); 
                c_1 = -1/(4*pi*delta^2);
                K_reg = @(x) c_0+c_1*x.^2;
                f_reg = @(x) 2*c_1*x;
            elseif reg == 2
                c_0 = -(-3+4*log(delta))/(8*pi);
                c_1 = -1/(2*pi*delta^2);
                c_2 = 1/(8*pi*delta^4);
                K_reg = @(x) c_0+c_1*x.^2+c_2*x.^4;
                f_reg = @(x) 2*c_1*x+4*c_2*x.^3;
            elseif reg == Inf
                K_reg = @(x) 0;
                f_reg = @(x) 0;
                c_0 = 0;
            end

            K = @(x) fillmissing(K_reg(x).*(abs(x)<delta)+K(x).*(abs(x)>=delta),'constant',c_0);
            f = @(x) fillmissing(f_reg(x).*(abs(x)<delta)+f(x).*(abs(x)>=delta),'constant',0);

            if and(p>0,mod(p,2) == 0)

                K = @(x) K(x) + x.^p/p;
                f = @(x) f(x) + x.^(p-1);

            elseif and(p>0,mod(p,2) ~= 0)

                K = @(x) K(x) + abs(x).^p/p;
                f = @(x) f(x) + sign(x).*abs(x).^(p-1);
            end

        else
      if and(p>0,mod(p,2) == 0)

                K = @(x) K(x) + x.^p/p;
                f = @(x) f(x) + x.^(p-1);

            elseif and(p>0,mod(p,2) ~= 0)

                K = @(x) K(x) + abs(x).^p/p;
                f = @(x) f(x) + sign(x).*abs(x).^(p-1);
            end
        end
    end
end
