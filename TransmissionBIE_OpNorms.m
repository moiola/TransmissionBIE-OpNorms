% 9.2021   TransmissionBIE_OpNorms.m    R.Hiptmair, A.Moiola, E.A.Spence
%
% Helmholtz transmission problems on unit disc and sphere via Fourier
% diagonalisation.
% Plot (against the wavenumber k) the norms of (i) the solution operator
% and (ii) the inverses of the 1st/2nd-kind classical/augmented BIOs.

function TransmissionBIE_OpNorms( Dim, ni, no )
tic; close all;

if nargin < 3
    Dim = 2;            % Space dimension, either 2 or 3.
    ni = 1;             % PDE coefficient in bounded domain (1 or 3 in paper example)
    no = 3;             % PDE coefficient in unbounded domain (3 or 1 in paper example)
end

kmax = 15;          % Maximal wavenumber k considered
numk = 1000;        % Number of values of k used (numk=10000 for the plots in the paper, reduce for quicker computations)
Lmax = 100;         % Maximal Fourier index considered (Lmax=100 for the plots in the paper)
SaveFig = 0;        % 1=Save figures plotted, 0 don't save

k = linspace(0.1,kmax,numk);
fprintf('\nPlot operator norms for %d-dimensional scattering problem, with ni=%g, no=%g, up to L=%d and k=%d, with %d values of k.\n', Dim,ni,no,Lmax,kmax,numk)
if ni > no, disp('ni>no: expect physical quasi-resonances.'),
else,       disp('ni<no: all quasi-resonances are spurious.'), end

%% Precompute Bessel/Hankel function values and BIO eigenvalues
if Dim == 2
    ell_vec = (-Lmax:Lmax)';        % Fourier indices
    Lshift = Lmax + 1;              % Shift to read correct Fourier indices from ell_vec
    [kk,ell] = meshgrid(k,ell_vec);
    ki = kk*sqrt(ni);
    ko = kk*sqrt(no);
    
    % Evaluations of Bessel J_l(kr) H^1_l(kr) and their derivatives at r=1
    % All these variables have dimension (2Lmax+1) x numk
    Ji = besselj(ell, ki);
    Jo = besselj(ell, ko);
    Hi = besselh(ell, 1, ki);
    Ho = besselh(ell, 1, ko);
    Ai = (ki/2) .* ( besselj(ell-1,ki) - besselj(ell+1,ki) );  % derivative of J_l(kr)
    Ao = (ko/2) .* ( besselj(ell-1,ko) - besselj(ell+1,ko) );
    Bi = (ki/2) .* ( besselh(ell-1,1,ki) - besselh(ell+1,1,ki) ); % derivative of H_l^1(kr)
    Bo = (ko/2) .* ( besselh(ell-1,1,ko) - besselh(ell+1,1,ko) );
    
    % BIO diagonalisation - from [Amini, Maines] Thm 2:
    Vi = 1i* pi/2 * ( Ji.*Hi );
    Vo = 1i* pi/2 * ( Jo.*Ho );
    Ki = 1i* pi/2 * ( Ji.*Bi ) + 0.5;   % Factor k is already included in derivative Bi
    Ko = 1i* pi/2 * ( Jo.*Bo ) + 0.5;
    Wi = -1i* pi/2 * ( Ai.*Bi );        % Our sign choice is opposite to [Amini, Maines]
    Wo = -1i* pi/2 * ( Ao.*Bo );
    
    SobolevWeight = sqrt(kk.^2+ell.^2);    % Factor in k-weighted norms
    Shape = 'circle';
    
elseif Dim == 3
    ell_vec = (0:Lmax)';
    Lshift = 1;
    [kk,ell] = meshgrid(k,ell_vec);
    ki = kk*sqrt(ni);
    ko = kk*sqrt(no);
    
    % Spherical Bessel functions 
    SphBessel = @(l,x) sqrt((pi/2)./x) .* besselj(l+0.5,x);
    SphHankel = @(l,x) sqrt((pi/2)./x) .* besselh(l+0.5,1,x);
    Ji = SphBessel(ell,ki);
    Jo = SphBessel(ell,ko);
    Hi = SphHankel(ell,ki);
    Ho = SphHankel(ell,ko);
    % Derivatives of  j_l(kr)  and h^1_l(kr)  at r=1 --- including factor k
    % using formula   f_n'(x)=f_{n-1}(x)-(n+1)/x f_n(x)   DLMF-10.51.2
    Ai = SphBessel(ell-1,ki) .* ki - (ell+1) .* Ji ;
    Ao = SphBessel(ell-1,ko) .* ko - (ell+1) .* Jo ;
    Bi = SphHankel(ell-1,ki) .* ki - (ell+1) .* Hi ;
    Bo = SphHankel(ell-1,ko) .* ko - (ell+1) .* Ho ;
   
    % BIOs following Table 1 in [Vico, Greengard, Gimbutas]:
    Vi = 1i * ki .* Ji .* Hi;    
    Vo = 1i * ko .* Jo .* Ho;
    Ki = 1i * 0.5 * ki .* (Ji.*Bi + Ai.*Hi);    
    Ko = 1i * 0.5 * ko .* (Jo.*Bo + Ao.*Ho);
    Wi = -1i * ki .* Ai .* Bi;      % Sign opposite to VGG
    Wo = -1i * ko .* Ao .* Bo;
    
    SobolevWeight = sqrt(kk.^2+ell.*(ell+1));   
    Shape = 'sphere';    
else
    error('Dimension must be 2 or 3')
end
Vi = Vi .* SobolevWeight;
Vo = Vo .* SobolevWeight;
Wi = Wi ./ SobolevWeight;
Wo = Wo ./ SobolevWeight;

%% Compute operator norms and check which Fourier indices give maximal values
NormAI = zeros(numk,1);
NormAII = zeros(numk,1);
NormSio = zeros(numk,1);
LargestL_AI = 0;
LargestL_AII = 0;
LargestL_Sio = 0;
for i_k = 1:numk
    NormAI_k = zeros(length(ell_vec),1);
    NormAII_k = zeros(length(ell_vec),1);
    NormSio_k = zeros(length(ell_vec),1);
    for i_l = 1:length(ell_vec)
    % Write operators for l-th Fourier mode as 2x2 matrices & compute norms
        AI = [-Ki(i_l,i_k)-Ko(i_l,i_k), Vi(i_l,i_k)+Vo(i_l,i_k); Wi(i_l,i_k)+Wo(i_l,i_k),Ki(i_l,i_k)+Ko(i_l,i_k)];
        AII = [1+Ki(i_l,i_k)-Ko(i_l,i_k), -Vi(i_l,i_k)+Vo(i_l,i_k); -Wi(i_l,i_k)+Wo(i_l,i_k),1-Ki(i_l,i_k)+Ko(i_l,i_k)];
        Sio = [-Ji(i_l,i_k)*Bo(i_l,i_k), Ji(i_l,i_k)*Ho(i_l,i_k).*SobolevWeight(i_l);...
               -Ai(i_l,i_k)*Bo(i_l,i_k)./SobolevWeight(i_l), Ai(i_l,i_k)*Ho(i_l,i_k)] / (Ho(i_l,i_k)*Ai(i_l,i_k)-Ji(i_l,i_k)*Bo(i_l,i_k));
        NormAI_k(i_l) = norm(inv(AI));
        NormAII_k(i_l) = norm(inv(AII));
        NormSio_k(i_l) = norm(Sio);
    end
    % For each k compute max norm over Fourier indices and store largest l
    [NormAI(i_k),Loc_AI] = max(NormAI_k);
    [NormAII(i_k),Loc_AII] = max(NormAII_k);
    [NormSio(i_k),Loc_Sio] = max(NormSio_k);
    LargestL_AI = max( LargestL_AI, abs(Loc_AI - Lshift));
    LargestL_AII = max( LargestL_AII, abs(Loc_AII - Lshift));
    LargestL_Sio = max( LargestL_Sio, abs(Loc_Sio - Lshift));    
end

%% Plot operator norms for classical (non-augmented) formulations
figure('position',[50,30,1100,750])
semilogy(k,NormAI, k,NormAII, 'linewidth',2)
hold on, grid on
semilogy(k,NormSio, 'linewidth',4);
xlabel('Wavenumber $k$','interpreter','latex')
ylabel('Operator norm', 'interpreter','latex')
title(['Unit ',Shape,' transmission problem, classical formulations: $n_i=',num2str(ni),'$, $n_o=',num2str(no),'$'],  'interpreter','latex')
legend('$\|A_I^{-1}\|$','$\|A_{II}^{-1}\|$','$\|S_{io}\|$', 'interpreter','latex', 'location','northwest')
addToolbarExplorationButtons(gcf); 
set(gca,'FontSize',22)
if SaveFig
    print('-dpng',['TransmissionBIE-OpNorms-',num2str(Dim),'d-ni',num2str(ni),'-no',num2str(no),'_k.png'])
end

%% Repeat for augmented formulation
if ni<no
    NormAIaug = zeros(numk,1);
    NormAIIaug = zeros(numk,1);
    LargestL_AIaug = 0;
    LargestL_AIIaug = 0;
    for i_k = 1:numk
        NormAIaug_k = zeros(length(ell_vec),1);
        NormAIIaug_k = zeros(length(ell_vec),1);
        for i_l = 1:length(ell_vec)
            AI = [-Ki(i_l,i_k)-Ko(i_l,i_k), Vi(i_l,i_k)+Vo(i_l,i_k); Wi(i_l,i_k)+Wo(i_l,i_k),Ki(i_l,i_k)+Ko(i_l,i_k)];
            AII = [1+Ki(i_l,i_k)-Ko(i_l,i_k), -Vi(i_l,i_k)+Vo(i_l,i_k); -Wi(i_l,i_k)+Wo(i_l,i_k),1-Ki(i_l,i_k)+Ko(i_l,i_k)];        
            AIaug = [AI; [.5+Ki(i_l,i_k), -Vi(i_l,i_k); -Wi(i_l,i_k), .5-Ki(i_l,i_k)]];
            AIIaug = [AII; [.5+Ki(i_l,i_k), -Vi(i_l,i_k); -Wi(i_l,i_k), .5-Ki(i_l,i_k)]];
            NormAIaug_k(i_l) = 1/min(svd(AIaug));       % Pseudo-inverse norm computed from SVD
            NormAIIaug_k(i_l) = 1/min(svd(AIIaug));
        end
        [NormAIaug(i_k),Loc_AIaug] = max(NormAIaug_k);
        [NormAIIaug(i_k),Loc_AIIaug] = max(NormAIIaug_k);
        LargestL_AIaug = max( LargestL_AIaug, abs(Loc_AIaug - Lshift));
        LargestL_AIIaug = max( LargestL_AIIaug, abs(Loc_AIIaug - Lshift));
    end
    
    figure('position',[300,30,1200,750])
    semilogy(k,NormAIaug, k,NormAIIaug, 'linewidth',3)
    hold on, grid on
    semilogy(k,NormSio, 'linewidth',4);
    xlabel('Wavenumber $k$','interpreter','latex')
    ylabel('Operator norm', 'interpreter','latex')
    title(['Unit ',Shape,' transmission problem, augmented formulations: $n_i=',num2str(ni),'$, $n_o=',num2str(no),'$'],  'interpreter','latex')
    legend('$\|\widetilde A_I^\dagger\|$','$\|\widetilde A_{II}^\dagger\|$','$\|S_{io}\|$', 'interpreter','latex', 'location','northwest')
    addToolbarExplorationButtons(gcf);
    set(gca,'FontSize',22)
    if SaveFig
        print('-dpng',['TransmissionBIE-OpNorms-Augmented-',num2str(Dim),'d-ni',num2str(ni),'-no',num2str(no),'_k.png'])
    end
    disp('The largest Fourier indices maximising the norms for some k are')
    fprintf(' A_I: %d\n A_II: %d\n S_io: %d\n A_I-aug: %d\n A_II-aug: %d\n',LargestL_AI,LargestL_AII, LargestL_Sio, LargestL_AIaug, LargestL_AIIaug)
else    
    disp('The largest Fourier indices maximising the norms for some k are')
    fprintf(' A_I: %d\n A_II: %d\n S_io: %d\n',LargestL_AI,LargestL_AII, LargestL_Sio)
end
toc
end