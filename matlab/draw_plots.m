function out = draw_plots(HVec, n_vec, V1, D, xt_h, xt, x_ph, OptTime)
% Draw the plots associated with REED_SDP.m

  figure(1), clf % PLOTTING FUNCTIONS
  subplot(2,2,1), hold on, set(gca,'fontsize',14)
  P=pcolor(HVec,n_vec,V1); set(P,'edgecolor','none'), colorbar, ylim([0 75])
  title('Reward in 1st timestep','fontsize',14)
  xlabel('Alternative harvest rates','fontsize',14)
  ylabel('Population','fontsize',14), axis tight
  subplot(2,2,2), cla, hold on, xlim([1 OptTime]), set(gca,'fontsize',14)
  ZZ = HVec(D);
  ZZ(ZZ==0)=nan;
  P=pcolor([OptTime:-1:1],n_vec,ZZ(:,end:-1:1));
  set(P,'edgecolor','none'), colorbar, ylim([0 75])

  title('Optimal decision (no scrap value)','fontsize',14)
  xlabel('Management timeline (years)','fontsize',14)
  ylabel('Population','fontsize',14), axis tight
  subplot(2,2,[3,4]), hold on, set(gca,'fontsize',14)

  ReedThreshold = n_vec(sum(D(:,1)==1));
  plot([0 OptTime],[ReedThreshold ReedThreshold],'r:')
  MixX = [1:length(xt);
  1:length(xt)];
  MixX = MixX(1:end-1);
  Mix = [xt;xt_h];
  Mix = Mix(1:end-1);
  plot(MixX,Mix,'g--','linewidth',2)
  plot([1:OptTime],xt_h,'linewidth',2)
  plot([1:OptTime],[x_ph xt_h(end)],'m--','linewidth',2), ylim([0 1.1*max([xt_h,xt])])
  xlabel('Management timeline (years)','fontsize',16)
  ylabel('Population','fontsize',16), axis tight
  legend('Reed''s S','N(t) without harvest','N(t)','Post-harvest population',2);
  %set(L,'fontsize',12)
  out = 0;
end

