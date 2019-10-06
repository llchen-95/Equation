% Newton-CG method for computing time-periodic and
% space-localized solutions in the CQGL equation:
% At-gamma*Axx+beta*|A|^2A+delta*|A|^4A-chi*A=0.
% Here, z represents scaled time tau, and A=U*exp(i*mu*t).
clc; clear all; close all;
load U0_fig3.mat;          % this data contains i.c. U(x, z)
Lx=100; Nx=512; Lz=2*pi; Nz=32; errormax=1e-10; errorCG=1e-4;
dx=Lx/Nx; x=-Lx/2:dx:Lx/2-dx; kx=[0:Nx/2-1  -Nx/2:-1]*2*pi/Lx;
dz=Lz/Nz; z=0:dz:Lz-dz;       kz=[0:Nz/2-1  -Nz/2:-1]*2*pi/Lz;
[X,Z]=meshgrid(x,z);  [KX,KZ]=meshgrid(kx,kz);  KX2=KX.*KX;

gamma=0.9-1.1i; beta=-3-i; delta=2.75-i; chi=-0.1;
gamma1=real(gamma); gamma2=imag(gamma);
beta1=real(beta);   beta2=imag(beta);
delta1=real(delta); delta2=imag(delta);

tic;
nnt=0;                       % nnt: # of Newton steps
ncg=0;                       % ncg: # of CG iterations
while 1                      % Newton-CG iterations start
  nnt=nnt+1;
  u=real(U); v=imag(U); U2=u.*u+v.*v; U4=U2.*U2;
  G=gamma*ifft2(-KX2.*fft2(U))-(beta*U2+delta*U4-chi).*U;
  Ut=ifft2(i*KZ.*fft2(U)); ut=real(Ut); vt=imag(Ut);
  produv=2*sum(sum(u.*v)); produtvt=2*sum(sum(ut.*vt));
  mu   = sum(sum(v.*imag(G)-u.*real(G)))/produv;
  omega= sum(sum(ut.*imag(G)+vt.*real(G)))/produtvt;
  L0U=omega*Ut+i*mu*U-G;
  Uerror(nnt)=max(max(abs(L0U))); Uerror(nnt)
  numcg(nnt)=ncg; time(nnt)= toc;
  if Uerror(nnt) < errormax
      break
  end

  betaU1=beta1*u-beta2*v;    betaU2=beta1*v+beta2*u;
  deltaU1=delta1*u-delta2*v; deltaU2=delta1*v+delta2*u;
  G11=chi-beta1*U2-betaU1*2.*u-delta1*U4-deltaU1*4.*u.*U2;
  G12=   +beta2*U2-betaU1*2.*v+delta2*U4-deltaU1*4.*v.*U2;
  G21=   -beta2*U2-betaU2*2.*u-delta2*U4-deltaU2*4.*u.*U2;
  G22=chi-beta1*U2-betaU2*2.*v-delta1*U4-deltaU2*4.*v.*U2;
  Dxx=@(F)   ifft2(-KX2.*fft2(F));
  Dtxx=@(F)  ifft2(( omega*i*KZ+gamma1*KX2).*fft2(F));
  DtxxA=@(F) ifft2((-omega*i*KZ+gamma1*KX2).*fft2(F));

  P=@(F)  Dtxx(real(F))-G11.*real(F)-(mu+G12).*imag(F)  ...
          +gamma2*Dxx(imag(F)) ...
     +i*( (mu-G21).*real(F)-gamma2*Dxx(real(F))   ...
          +Dtxx(imag(F))-G22.*imag(F) );
  PA=@(F) DtxxA(real(F))-G11.*real(F)+(mu-G21).*imag(F) ...
          -gamma2*Dxx(imag(F)) ...
     +i*( -(mu+G12).*real(F)+gamma2*Dxx(real(F))  ...
          +DtxxA(imag(F))-G22.*imag(F) );

  L1= @(F) P(F)-sum(sum(imag(Ut.*P(F))))/produtvt*Ut ...
              +sum(sum(real(U.*P(F))))/produv*i*U;
  L1A=@(F) PA(F)  ...
     -sum(sum(real(conj(F).*Ut)))/produtvt*PA(vt+i*ut) ...
     -sum(sum(imag(conj(F).*U)))/produv*PA(u-i*v);

  c=8;                   % Here fftM is the preconditioner
  fftM=omega^2*KZ.*KZ+abs(gamma)^2*KX2.*KX2+c;
  dU=0*Z;                            % CG iterations start
  R=-L1A(L0U);
  MinvR=ifft2(fft2(R)./fftM);
  R2=sum(sum(real(conj(R).*MinvR))); R20=R2;
  D=MinvR;
  while (R2 > R20*errorCG^2)
      L2D=L1A(L1(D));
      a=R2/sum(sum(real(conj(D).*L2D)));
      dU=dU+a*D;
      R=R-a*L2D;
      MinvR=ifft2(fft2(R)./fftM);
      R2old=R2;
      R2=sum(sum(real(conj(R).*MinvR)));
      b=R2/R2old;
      D=MinvR+b*D;
      ncg=ncg+1;
  end                                  % CG iterations end
  U=U+dU;
end                             % Newton-CG iterations end

% plotting of numerical results
subplot(221); imagesc(x, [z z+Lz]/omega, abs([U;U]));
axis xy; colorbar; xlabel('x'); ylabel('t'); title('(a)');
subplot(222); imagesc(x, [z z+Lz]/omega, angle([U;U]));
axis xy; colorbar; xlabel('x'); ylabel('t'); title('(b)');
subplot(223); semilogy(numcg, Uerror, numcg, Uerror, 'o');
xlabel('number of CG iterations'); ylabel('solution error');
title('(c)');
subplot(224); semilogy(time/60, Uerror, time/60, Uerror, 'o');
xlabel('time (minutes)'); ylabel('solution error');
title('(d)');
format long; period=2*pi/omega
mu