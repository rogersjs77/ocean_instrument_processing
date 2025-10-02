function [F] = get_wave_friction(A,B,dx)

rho = 1022;
g = 9.81;
F.dx = dx;

% get C, Cg for time and freqency
for i =1:size(A.DIR,1)
[A.C(i,:),A.Cg(i,:)] = get_c_cg(2*pi*A.f,A.h(i));
[B.C(i,:),B.Cg(i,:)] = get_c_cg(2*pi*B.f,B.h(i));

% get k
A.k(i,:) = get_wavenumber(2*pi*A.f,A.h(i));
B.k(i,:) = get_wavenumber(2*pi*B.f,B.h(i));
end

% compute spectral quantities
for i =1:size(A.DIR,1)
%  spectral energy flux
   A.Flux(i,:)=rho*g*(A.Cg(i,:).*cosd(A.DIR(i,:)).*A.SSE(i,:)); 
   B.Flux(i,:)=rho*g*(B.Cg(i,:).*cosd(B.DIR(i,:)).*B.SSE(i,:));   
       
%  spectral urms squared
   A.Urms_sqrd(i,:) = 4*pi^2*A.f.^2.*A.SSE(i,:)...
       ./(sinh(A.k(i,:).*A.h(i))).^2;
   B.Urms_sqrd(i,:) = 4*pi^2*B.f.^2.*B.SSE(i,:)...
       ./(sinh(B.k(i,:).*B.h(i))).^2;
   
%  spectral Ab
   A.Ab_sqrd(i,:) = 2*(A.SSE(i,:)./(sinh(A.k(i,:).*A.h(i))).^2);
   B.Ab_sqrd(i,:) = 2*(B.SSE(i,:)./(sinh(B.k(i,:).*B.h(i))).^2);
   
end

F.df = A.f(2)-A.f(1);
indx = A.f < 1/5 & A.f > 1/25;

% spectral sum quantities, over indx
F.Flux(1,:) = nansum(A.Flux(:,indx).*F.df,2);
F.Flux(2,:) = nansum(B.Flux(:,indx).*F.df,2);

F.Hs(1,:) = 4*sqrt(nansum(A.SSE(:,indx)*F.df,2));
F.Hs(2,:) = 4*sqrt(nansum(B.SSE(:,indx)*F.df,2));

F.Urms(1,:) = sqrt(nansum(A.Urms_sqrd(:,indx).*F.df,2));
F.Urms(2,:) = sqrt(nansum(B.Urms_sqrd(:,indx).*F.df,2));

F.Ab(1,:) = sqrt(nansum(A.Ab_sqrd(:,indx).*F.df,2));
F.Ab(2,:) = sqrt(nansum(B.Ab_sqrd(:,indx).*F.df,2));

% average quantities
F.Urms_cubed_avg_spc = 0.5*(A.Urms_sqrd.^1.5+B.Urms_sqrd.^1.5);
F.Urms_cubed_avg = 0.5*(F.Urms(1,:).^3+F.Urms(2,:).^3);
F.Ab_avg = 0.5*(F.Ab(1,:)+F.Ab(2,:));

% calculated bottom dissipation
F.dissip_calc = (F.Flux(1,:)-F.Flux(2,:))/F.dx;

% compute fw in time
F.fw_calc = F.dissip_calc./(0.6*rho*F.Urms_cubed_avg);
F.fw_calc(F.dissip_calc<1) = nan;

% find constant fw for lowest RMSE
F.fw=0.1;
rmse1=inf;
drmse=-1;
i=0;
while drmse<0 && i<1E5
    F.dissip_model_const = 0.6*F.fw*rho*F.Urms_cubed_avg;
    [F.r2_const, rmse2]=rsquare(F.dissip_model_const,F.dissip_calc);
    F.fw=1.001*F.fw;
    drmse = rmse2-rmse1;
    rmse1=rmse2;
    i=i+1;
end
F.rmse_const = rmse1;
[F.skill_const,~] = skill(F.dissip_model_const,F.dissip_calc);
clear drmse rmse1 rmse2 i

% find kN for lowest RMSE using Nielsen 92 model
F.kN = 0.1;
rmse1=1E4;
drmse = -1;
i=0;
while drmse<0 && i<1E5
    % nielsen 92 method
%     a1=5.5;
%     a2=-0.20;
%     a3 = -6.3;
    % swart 74(Jonnson 66)
    a1=5.213;
    a2=-0.194;
    a3=-5.977;
    fw_n92=exp(a1*(F.Ab_avg/F.kN).^a2+a3);
%     fw_n92(Ab_avg/kN<0.172) = 4.58;

    F.dissip_model_n92 = 0.6*fw_n92.*rho.*F.Urms_cubed_avg;
    [F.r2_variable, rmse2]=rsquare(F.dissip_model_n92,F.dissip_calc);
    F.kN=1.001*F.kN;
    drmse = rmse2-rmse1;
    rmse1=rmse2;
    i=i+1;
end
F.rmse_variable = rmse1;
[F.skill_variable,~] = skill(F.dissip_model_n92,F.dissip_calc);

clear drmse rmse1 rmse2 i


%% calculate fw spectrally
F.fw_spc_calc = (A.Flux-B.Flux)./ (F.dx*0.6*rho*F.Urms_cubed_avg_spc);
% F.fw_spc_calc(F.dissip_calc<5E-4)=nan;
F.fw_spc_calc_avg = nanmean(F.fw_spc_calc,1);
F.fw_spc_calc_avg(F.dissip_calc<1)=nan;

%% export variables
F.A = A;
F.B = B;

end
