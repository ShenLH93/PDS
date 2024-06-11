function [I, R ,t,p] = monoPDS(Imax ,Imin, alpha, beta, lambda, vareps, r, r0, K, debug)

%% Tradeoff parameters
if (~exist('alpha','var'))	% alpha -- parameter for shape
    alpha = 0.001;
end
if (~exist('beta','var'))	% beta -- parameter for texture
    beta = 0.0001;
end
if (~exist('lambda','var'))	% lambda -- parameter for illumination
    lambda = 0.25;
end
if (~exist('vareps','var')) % vareps -- stopping parameter
    vareps = 0.001;
end
if (~exist('K','var'))      % K -- maximum iterations
    K = 5;
end
if (~exist('r','var'))      % r -- the size of Omega in Eq.(3)
    r = 5;
end
if (~exist('r0','var'))     % r0 -- the size of Omega in Eq.(7)
    r0 = 3;
end
if (~exist('debug','var'))  % debug -- set debug/release
    debug = true;
end


r = (r-1)/2;
r0 = (r0-1)/2;
eps=0.0001;


%% Initialize the original scene intensity S,
%% the intensity of the scene at infinity A∞
%% and degree of liner polarization p.
Iminus = Imax-Imin;
Isum   = Imax+Imin;
S = Isum;
A = convMax(single(S),r0);
A = guidedfilter(S, A, 20, eps);
p = Iminus./Isum;

%% initialize backward scattering Bg and direct transmission D_tmp
e = 1.13;  %  Yoav Y. Schechner 'Recovery of Underwater Visibility and Structure by Polarization Analysis'
Bg = Iminus./(p*e);
D_tmp = Isum-Bg;
D_tmp = max(D_tmp,0);

%% initialize transmission map t
Bmax = min(max(max(Isum)),0.9);
t = 1-Isum.*Iminus./Bmax./1.13;
t = max(t,0.01);
t = convMax(single(t),r0);
t= double(t);
%% initialize the ambient illumination L
L = D_tmp./t;

%% initialize Illumination I_0
I =L;
%% initialize reflectance R_0
R=ones(size(S));

%% Optimization Process 
if debug == true
    fprintf('-- Stop iteration until eplison < %02f or K > %d\n', vareps, K);
end

for iter = 1:K
    preI = I;
    preR = R;
    pret = t;

    %% algorithm for optimizing Ik 
    I=(S-(1-t).*Bmax)./R./t;
    Ix = diff(I,1,2); Ix = padarray(Ix, [0 1], 'post');
    Iy = diff(I,1,1); Iy = padarray(Iy, [1 0], 'post');
    avgIx=convBox( single(Ix), r);
    avgIy=convBox( single(Iy), r);
    ux = max(abs(avgIx.*Ix),eps).^(-1);% ux in Eq.(16)
    uy = max(abs(avgIy.*Iy),eps).^(-1);% uy in Eq.(16)
    ux(:,end) = 0;
    uy(end,:) = 0;

    Index1a = S-(1-t).*Bmax;
    Index1b = R.*t;

    I = solveLinearSystem(Index1a, Index1b, ux, uy, alpha, A, lambda);  % Eq.(21)
    eplisonI = norm(I-preI, 'fro')/norm(preI, 'fro');% iterative error of I

    %% algorithm for optimizing Rk
    R=(S-(1-t).*Bmax)./I./t;
    Rx = diff(R,1,2); Rx = padarray(Rx, [0 1], 'post');
    Ry = diff(R,1,1); Ry = padarray(Ry, [1 0], 'post');
    vx = max(abs(Rx),eps).^(-1);  % vx in Eq.(16)
    vy = max(abs(Ry),eps).^(-1);  % vy in Eq.(16)
    vx(:,end) = 0;
    vy(end,:) = 0;
    Index2a = S-(1-t).*Bmax;
    Index2b  = I.*t;
    R = solveLinearSystem(Index2a, Index2b, vx, vy, beta); % Eq.(23)
    eplisonR = norm(R-preR, 'fro')/norm(preR, 'fro'); % iterative error of R
    
    %% algorithm for optimizing Iminus
    Iminus = guidedfilter(R, Iminus, 3, eps);

    %% algorithm for optimizing tk 
    t = 1-Isum.*Iminus./Bmax./1.13;
    tx = diff(t,1,2); tx = padarray(tx, [0 1], 'post');
    ty = diff(t,1,1); ty = padarray(ty, [1 0], 'post');
    wx = max(abs(tx),eps).^(-1); 
    wy = max(abs(ty),eps).^(-1); 
    wx(:,end) = 0;
    wy(end,:) = 0;
    Index3a = 1-Isum.*Iminus./Bmax./1.13;
    Index3b = ones(size(I));
    t = solveLinearSystem(Index3a, Index3b, wx, wy, beta);

    eplisont = norm(t-pret, 'fro')/norm(pret, 'fro'); % iterative error of t
    
    if mean(mean(t))>0.95
        t= ones(size(I));
        eplisont = 1;
    end

    %% Calculate the optimized degree of polarization p
    p = Iminus./Isum;

    %% iteration until convergence
    if debug == true
        fprintf('Iter #%d : eplisonI = %f; eplisonR = %f; eplisont = %f\n', ...
            iter, eplisonI, eplisonR, eplisont);
    end
    if(eplisonI<vareps||eplisonR<vareps||eplisont<vareps*0.1*vareps)
        break;
    end

end
end

