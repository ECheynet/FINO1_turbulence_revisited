function [Coeff,S5] = fit_and_plot_spectra(Su,eta,target,minU,Utot,Wind,z,u_star,component,guess1)
%FIT_AND_PLOT_SPECTRA Fits turbulence spectra based on stability classes and plots results.
%
%   [Coeff, S5] = FIT_AND_PLOT_SPECTRA(Su, eta, target, minU, Utot, Wind, z, u_star, component, guess1)
%
%   Inputs:
%     Su        : Spectral density (size: [3 x N x f]) for the velocity component or co-spectrum.
%     eta       : Stability parameter z/L for each sample and height.
%     target    : Vector of bin edges for stability classification.
%     minU      : Minimum wind speed used for frequency normalization.
%     Utot      : Mean wind speed matrix [3 x N].
%     Wind      : Struct with fields: Wind.f (frequencies), Wind.time, etc.
%     z         : Measurement heights [3 x 1] in meters.
%     u_star    : Friction velocity array [3 x N].
%     component : Character ('u','v','w','uw') specifying the velocity component.
%     guess1    : Initial guess vector for non-linear least squares fitting (4 parameters).
%
%   Outputs:
%     Coeff     : Struct with fitted model parameters per stability bin and height.
%     S5        : Struct with median, 10th/90th percentiles, and Mann-scaled spectra for near-neutral cases.
%
% Loop over stability bins and heights
%   - Filter data by stability bin and height
%   - Normalize frequency and spectra
%   - Compute medians and quantiles
%   - Fit spectral models (non-linear least squares)
%   - Choose best model based on residuals or bounds
%   - Store results and plot fits vs measurements
%
% Output Coeff contains model parameters for 3 possible models (C1, C2, C3)
%
% The plotting is done using tiledlayout to visualize spectra for all stability bins
% Reference curves (Kaimal models) are included for comparison under neutral conditions
%
% This function reproduces core results from:
% Cheynet et al. (2018) "Velocity spectra and coherence estimates in the marine ABL"
% BLM, 169, 429–460. DOI: 10.1007/s10546-018-0382-2
%
% --------------------
% Author: Etienne Cheynet
% Contact: etienne.cheynet@uib.no
% --------------------

%% Inputparseer
if strcmpi(component,'uw'),
    myFun1 = @(a,k) a(1)*k(2,:)./(1+a(2).*k(1,:)).^(7/3)+a(3)*k(2,:)./(1+a(4).*k(1,:).^(7/3));
    myFun2 = @(a,k) a(1)*k(2,:)./(1+a(2).*k(1,:)).^(7/3)+a(3)*k(1,:)./(1+a(4).*k(1,:).^(7/3))+a(5).*k(2,:).^(-2);
    myFun3 = @(a,k) a(1).*k(2,:).^(-4/3)+a(2)*k(1,:)./(1+a(3).*k(1,:).^(7/3))+a(4).*k(2,:).^(-2);

    myFun1a = @(a,k) a(1)*k./(1+a(2).*k).^(7/3)+a(3)*k./(1+a(4).*k.^(7/3));
    myFun2a = @(a,k) a(1)*k./(1+a(2).*k).^(7/3)+a(3)*k./(1+a(4).*k.^(7/3))+a(5).*k.^(-2);
    myFun3a = @(a,k) a(1).*k.^(-4/3)+a(2)*k./(1+a(3).*k.^(7/3))+a(4).*k.^(-2);
elseif strcmpi(component,'u')||strcmpi(component,'v')||strcmpi(component,'w')
    myFun1 = @(a,k) a(1)*k(2,:)./(1+a(2).*k(1,:)).^(5./3)+a(3)*k(2,:)./(1+a(4).*k(1,:).^(5/3));
    myFun2 = @(a,k) a(1)*k(2,:)./(1+a(2).*k(1,:)).^(5/3)+a(3)*k(1,:)./(1+a(4).*k(1,:).^(5/3))+a(5).*k(1,:).^(-2);
    myFun3 = @(a,k) a(1).*k(2,:).^(-2/3)+a(2)*k(1,:)./(1+a(3).*k(1,:).^(5/3))+a(4).*k(1,:).^(-2);

    myFun1a = @(a,k) a(1)*k./(1+a(2).*k).^(5./3)+a(3)*k./(1+a(4).*k.^(5/3));
    myFun2a = @(a,k) a(1)*k./(1+a(2).*k).^(5/3)+a(3)*k./(1+a(4).*k.^(5/3))+a(5).*k.^(-2);
    myFun3a = @(a,k) a(1).*k.^(-2/3)+a(2)*k./(1+a(3).*k.^(5/3))+a(4).*k.^(-2);
end

Coeff.myFun1 = myFun1a;
Coeff.myFun2 = myFun2a;
Coeff.myFun3 = myFun3a;

if strcmpi(component,'u')||strcmpi(component,'v')
    guess2 = [10,4,5,1,0.01];
    guess3 = [4,5,1,0.1];
elseif strcmpi(component,'w')
    guess2 = [1,1,1,1,1e-12];
    guess3 = [1,1,1,1e-12];
elseif strcmpi(component,'uw')
    guess2 = [1,1,1,1,1e-6];
    guess3 = [10,1,1,1e-6];
end


figure('position',[100 500 1200 600])

options=optimset('TolX',1e-8,'TolFun',1e-8,'Display','off');
Nbin = 60;
Coeff.C1 = nan(numel(target)-1,3,numel(guess1));
Coeff.C2 = nan(numel(target)-1,3,numel(guess2));
Coeff.C3 = nan(numel(target)-1,3,numel(guess3));
Coeff.zL = nan(numel(target)-1,numel(z));

tiledlayout(3,3,"TileSpacing","tight")
for ii=1:numel(target)-1
    nexttile
    grid on
    N = zeros(1,3);
    for indZ = 1:3
        hold on;box on;
        indL = find(eta(indZ,:)>=target(ii) & eta(indZ,:)<target(ii+1));
        if numel(indL)>1,
            Coeff.zL(ii,indZ) = median(eta(indZ,indL),'omitnan');
            N(indZ) = numel(indL);

            clear n
            n(1,:) = logspace(log10(1/3600*z(indZ)/minU),log10(10*z(indZ)/nanmedian(Utot(indZ,:))),Nbin);
            n(2,:) = logspace(log10(1/3600*z(indZ)/minU),log10(10*z(indZ)/nanmedian(Utot(indZ,:))),Nbin);

            if ii>6,
                xlabel('$nz/\overline{u}$','interpreter','latex')
            end

            if ii ==1 || ii ==4 || ii ==7
                if strcmpi(component,'u'),
                    ylabel('$n S_u/u^2_*$','interpreter','latex')
                elseif strcmpi(component,'w'),
                    ylabel('$n S_w/u^2_*$','interpreter','latex')
                elseif strcmpi(component,'v'),
                    ylabel('$n S_v/u^2_*$','interpreter','latex')
                elseif strcmpi(component,'uw'),
                    ylabel('$n |Co_{uw}|/u^2_*$','interpreter','latex')
                end
            end



            X0 = bsxfun(@times,1./Utot(indZ,indL),Wind.f.*z(indZ))';
            X1 = bsxfun(@times,1./Utot(indZ,indL),Wind.f.*z(indZ))';


            S = squeeze(Su(indZ,indL,:));
            Stot = S.*(Wind.f*u_star(indZ,indL).^(-2))';

            clear X
            X(1,:) = median(X0,'omitnan');
            X(2,:) = median(X1,'omitnan');
            Stot = Stot(:,~isnan(X(1,:)));
            X = X(:,~isnan(X(1,:)));
            if strcmpi(component,'w')|strcmpi(component,'uw'),
                indX = find(X(1,:) <12 );
            else
                indX = find(X(1,:) <X(1,end-14));
            end


            if abs(median(eta(indZ,indL),'omitnan'))<0.1 || numel(target)==2;
                S5.Su50(indZ,:) = median(Stot,'omitnan');
                S5.Su50_Mann(indZ,:) =  median(2*pi*X0(:,~isnan(X(1,:))).*S(:,~isnan(X(1,:)))./z(indZ),'omitnan');
                S5.Su90(indZ,:) = quantile(Stot,0.9);
                S5.Su10(indZ,:) = quantile(Stot,0.1);
                S5.f(indZ,:) = X(1,:);
                S5.n = Wind.f;
            elseif ii == numel(target)-1 && abs(nanmedian(eta(indZ,indL)))>=0.1 && ~exist('S5','var'),
                S5 = nan;
            end

            lim1 = zeros(1,numel(guess1));
            lim2 = [3e3*ones(1,numel(guess1))];
            [myPara50,~] = lsqcurvefit(myFun1,guess1,X(:,indX),median(Stot(:,indX),'omitnan'),lim1,lim2,options);
            if (any(myPara50>=500) && strcmpi(component,'w'))
                guess1 = [0.1,0.1,2,5];
                [myPara50,~] = lsqcurvefit(myFun1,guess1,X(:,indX),median(Stot(:,indX),'omitnan'),lim1,lim2,options);
                if (any(myPara50>=500) && strcmpi(component,'w'))
                    guess1 = [0,0,1,1];
                    [myPara50,~] = lsqcurvefit(myFun1,guess1,X(:,indX),median(Stot(:,indX),'omitnan'),lim1,lim2,options);
                end

            end

            if (any(myPara50>=500) && strcmpi(component,'u'))||(median(eta(indZ,indL),'omitnan')>-0.1 && strcmpi(component,'v')),
                Coeff.C3(ii,indZ,:) = lsqcurvefit(myFun3,guess3,X(:,indX),median(Stot(:,indX),'omitnan'),[0,0,0,0],[101,101,101,10],options);
                flag = 3; % trigger model 3
            else
                Coeff.C1(ii,indZ,:)= myPara50;
                flag = 1;
            end


            myX = X(1,:);

            if indZ ==1,
                g(1) = plot(myX,median(Stot,'omitnan'),'b^','markerfacecolor','c','markersize',3);
                if flag == 1,
                    h(1) = plot(myX,myFun1(Coeff.C1(ii,indZ,:),X),'b');
                elseif flag == 2,
                    h(1) = plot(myX,myFun2(Coeff.C2(ii,indZ,:),X),'b');
                elseif flag == 3,
                    h(1) = plot(myX,myFun3(Coeff.C3(ii,indZ,:),X),'b');
                end

            elseif indZ ==2,
                g(2) = plot(myX,median(Stot,'omitnan'),'ro','markerfacecolor',[1,0.6,0.8],'markersize',3);

                if flag == 1,
                    h(2) = plot(myX,myFun1(Coeff.C1(ii,indZ,:),X),'color','r','linestyle','-');
                elseif flag == 2,
                    h(2) = plot(myX,myFun2(Coeff.C2(ii,indZ,:),X),'color','r','linestyle','-');
                elseif flag == 3,
                    h(2) = plot(myX,myFun3(Coeff.C3(ii,indZ,:),X),'color','r','linestyle','-');
                end

            elseif indZ ==3,
                g(3) = plot(myX,median(Stot,'omitnan'),'kd','markerfacecolor','k','markersize',3);

                if flag == 1,
                    h(3) = plot(myX,myFun1(Coeff.C1(ii,indZ,:),X),'color','k');
                elseif flag == 2,
                    h(3) = plot(myX,myFun2(Coeff.C2(ii,indZ,:),X),'color','k');
                elseif flag == 3,
                    h(3) = plot(myX,myFun3(Coeff.C3(ii,indZ,:),X),'color','k');
                end


                for pp=1:numel(h),            uistack(h(pp), 'top'); end
                myEta = (target(ii)+target(ii+1))/2;
                if strcmpi(component,'u'),
                    ylim([0.01,5])
                    %                     xlim([1.1e-3,35])
                    S_kaimal = 105.*n(1,:).*(1+33.*n(1,:)).^(-5/3);
                    if abs(myEta)<0.3,
                        legKaimal = plot(n(1,:),S_kaimal,'k--');
                    else
                        legKaimal = plot(n(1,:),nan(size(S_kaimal)),'k--');
                    end
                elseif strcmpi(component,'v'),
                    ylim([0.005,10])
                    %                     xlim([1.1e-3,35])
                    S_kaimal = 17.*n(1,:).*(1+9.5.*n(1,:)).^(-5/3);
                    if abs(myEta)<0.3,
                        legKaimal = plot(n(1,:),S_kaimal,'k--');
                    else
                        legKaimal = plot(n(1,:),nan(size(S_kaimal)),'k--');
                    end
                elseif strcmpi(component,'w'),
                    ylim([0.001,2])
                    %                     xlim([1.1e-3,40])
                    S_kaimal = 2.1.*n(1,:)./(1+5.3.*n(1,:).^(5/3));
                    if abs(myEta)<0.3,
                        legKaimal = plot(n(1,:),S_kaimal,'k--');
                    else
                        legKaimal = plot(n(1,:),nan(size(S_kaimal)),'k--');
                    end
                elseif strcmpi(component,'uw'),
                    ylim([0.001,2])
                    %                     xlim([1.1e-3,40])
                    S_kaimal = 14.*n(1,:).*(1+9.6.*n(1,:)).^(-2.4);
                    if abs(myEta)<0.3,
                        legKaimal = plot(n(1,:),S_kaimal,'k--');
                    else
                        legKaimal = plot(n(1,:),nan(size(S_kaimal)),'k--');
                    end
                end
                if abs(myEta)<0.3,            uistack(legKaimal, 'top'); end

            end

            if ii<7,        set(gca,'xticklabel',{[]}) ;    end

            if ~any(ii==[1,4,7]),   set(gca,'yticklabel',{[]});    end

            if ii==2 && indZ ==3,
                leg = legend([g,legKaimal,h],'Measured ($z = 80$ m)',...
                    'Measured ($z = 60$ m)',...
                    'Measured ($z = 40$ m)',...
                    'Kaimal',...
                    'Fitted ($z = 80$ m)',...
                    'Fitted ($z = 60$ m)','Fitted ($z = 40$ m)','location','best');
                set(leg,'interpreter','latex','orientation','horizontal','Position',[0.2 0.94 0.6102 0.0213]);
            end

        end
        set(gca,'xscale','log','yscale','log');

        if indZ ==3,
            titlestring = [num2str(target(ii)),' $\leq \zeta <$',num2str(target(ii+1))];
            fth = text(.10,0.30,titlestring,...
                'units','normalized',...
                'horizontalalignment','left',...
                'verticalalignment','bottom','fontsize',10);
            set(gcf,'CurrentAxes',gca,'name',titlestring);
            set(fth,'interpreter','latex')

            titlestring = ['$N $($z=80$ m) = ',num2str(N(1))];
            fth = text(.1,0.17,titlestring,...
                'units','normalized',...
                'horizontalalignment','left',...
                'verticalalignment','bottom','fontsize',10);
            set(gcf,'CurrentAxes',gca,'name',titlestring);
            set(fth,'interpreter','latex')

            titlestring = ['$N $($z=60$ m) = ',num2str(N(2))];
            fth = text(.1,0.10,titlestring,...
                'units','normalized',...
                'horizontalalignment','left',...
                'verticalalignment','bottom','fontsize',10);
            set(gcf,'CurrentAxes',gca,'name',titlestring);
            set(fth,'interpreter','latex')

            titlestring = ['$N $($z=40$ m) = ',num2str(N(3))];
            % Make a title:
            fth = text(.1,0.03,titlestring,...
                'units','normalized',...
                'horizontalalignment','left',...
                'verticalalignment','bottom','fontsize',10);
            set(gcf,'CurrentAxes',gca,'name',titlestring);
            set(fth,'interpreter','latex')
        end

        set(gca,'TickLabelInterpreter','latex')
    end
end


set(gca,'GridLineStyle','-','GridColor',[0.7,0.7,0.7])
set(gcf,'color','w')

end


