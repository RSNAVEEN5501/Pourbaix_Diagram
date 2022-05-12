clc
figure('Name','Draft Pourbaix Diagram','NumberTitle','off');
hold on
prompt0 = {'Lower Limit of X','Upper Limit of X','Lower Limit of Y','Upper Limit of Y'};
data0 = {'0','14','-3','3'};
p0 = inputdlg(prompt0,'Axis',[1,40],data0);
dp0 = str2double(p0);
axis([dp0(1) dp0(2) dp0(3) dp0(4)])

fprintf('Note 1: Enter Gibbs Free Energy values in terms of J/mol');
fprintf('\nNote 2: If Hydrides are given then replace it with Hydroxides');
fprintf('\nNote 3: If equation has only e- then it is a horizontal line');
fprintf('\nNote 4: If equation has only H+ ions then it is a vertical line');
fprintf('\nNote 5: If equation has both e- and H+ then it is a sloped line');

%Initital Values%
prompt1 = {'Number of Reactions','Temperature (K)','Concentration of Ion','G of Metal','G of Metal Ion 1','G of Metal Ion 2','G of Metal Oxide 1','G of Metal Oxide 2','G of Metal Oxide Ion','G of Metal Hydroxide','G of Alloy Oxide','G of Water','Partial Pressure of Gases in atm','No of Horizontal Line','No of Vertical Line','No of Sloped Line'};
data1 = {'5','373','1e-06','-10767','-409935','0','-2017975','0','-916882','0','0','-312520','1','1','2','2'};
p1 = inputdlg(prompt1,'Initial Values',[1,40],data1);
dp1 = str2double(p1);
noe = dp1(1); temp = dp1(2); Conc = dp1(3); g_m = dp1(4); g_mi1 = dp1(5); g_mi2 = dp1(6); 
g_mo1 = dp1(7); g_mo2 = dp1(8); g_moi = dp1(9); g_mho = dp1(10); g_m1m2o=dp1(11); g_w = dp1(12);
g_h=0; g_e=0; ppog = dp1(13); horil = dp1(14); vertil = dp1(15); slopl = dp1(16);
F = 96485; R = 8.314;
c1 = ((2.303*R*temp)/F);

%Storage Values
f1 = 1; f2 = 1; f3 = 1;
Ea = zeros(horil,1);
pH = zeros(vertil,1);
EwpH = zeros(slopl,15);

%Constant Water Stability%
X_axis = 0:1:14;
E_H2 = -(0.0591*X_axis);
E_O2 = 1.23-(0.0591*X_axis);
plot(X_axis,E_H2,'--.r',X_axis,E_O2,'--.r');
xlabel('Potential of Hydrogen (pH) ----->');
ylabel('Potential (V) ----->');

fprintf('\nNote 6: Enter your equations with electrons and protons on reacting side');
if(noe==0)
    fprintf('\nPlease Enter the Number of Equations!');
end

%Co-Efficient
if(noe>=1)
    for i1=1:noe
        prompt2 = {'Coeff of Reactant Metal Ion 1','Coeff of Reactant Metal Ion 2','Coeff of Reactant Metal Oxide 1','Coeff of Reactant Metal Oxide 2','Coeff of Reactant Metal Oxide Ion','Coeff of Reactant Metal Hydroxide','Coff of Reactant Alloy Oxide','Coeff of Reactant H+','Coeff of Reactant e-','Coeff of Product Metal','Coeff of Product Metal Ion 1','Coeff of Product Metal Ion 2','Coeff of Product Metal Oxide 1','Coeff of Product Metal Oxide 2', 'Coeff of Product Metal Oxide Ion','Coeff of Product Metal Hydroxide','Coff of Product Alloy Oxide', 'Coeff of Product Water'};
        p2 = inputdlg(prompt2, "Reactions", [1,50]);
        dp2 = str2double(p2);
        a=dp2(1);b=dp2(2);c=dp2(3);d=dp2(4);e=dp2(5);f=dp2(6);g=dp2(7);m=dp2(8);n=dp2(9);
        h=dp2(10);i=dp2(11);j=dp2(12);k=dp2(13);l=dp2(14);o=dp2(15);p=dp2(16);q=dp2(17);r=dp2(18);
        %Assigning%
        arr1 = [a b c d e f g m n];
        arr2 = [h i j k l o p q r];
        comarr1 = [g_mi1 g_mi2 g_mo1 g_mo2 g_moi g_mho g_m1m2o g_h g_e];
        comarr2 = [g_m g_mi1 g_mi2 g_mo1 g_mo2 g_moi g_mho g_m1m2o g_w];
        
        %Calculation of Gibbs Free Energy%
        greact=0;
        gpdt=0;
        for i2=1:length(arr1)
            if(arr1(i2)>=1)
                greact = greact + (arr1(i2)*comarr1(i2));
            end
        end
        for i3=1:length(arr2)
            if(arr2(i3)>=1)
                gpdt = gpdt + (arr2(i3)*comarr2(i3));
            end
        end
        totg = gpdt-greact;

        %Calculation of Concentration%
        brr1 = [a b e];
        brr2 = [i j o];
        logr=0;
        logp=0;
        for i4=1:length(brr1)
            logr = logr + (brr1(i4)*log10(Conc)); 
        end
        for i5=1:length(brr2)
            logp = logp + (brr2(i5)*log10(Conc)); 
        end
        totlog = logp-logr;

        %Partial Pressure Calculation%
        if(temp==298)
            tppog = 0;
        elseif(temp>298 || temp<298)
            crr1 = [0 0];
            crr2 = [r 0];
            ppogr = 0;
            ppogp = 0;
            for i6=1:length(crr1)
                ppogr = ppogr + (crr1(i6)*log10(ppog)); 
            end
            for i7=1:length(crr2)
                ppogp = ppogp + (crr2(i7)*log10(ppog)); 
            end
            tppog = ppogp-ppogr;
        end

        %Conditions%
        if((m==0)&&(n>=1))
            enot1 = totg/(-n*F);
            ctemp1 = c1/n;
            tote1 = enot1-(ctemp1*totlog)-(ctemp1*tppog);
            Ea(f1,:) = tote1;
            f1=f1+1;
            yline(tote1,'-b');
        elseif((m>=1)&&(n==0))
            Keq = exp((-totg)/(R*temp));
            Klog = log10(Keq);
            pH1 = (Klog-tppog-totlog)/m;
            pH(f2,:) = pH1;
            f2=f2+1;
            xline(pH1,'-k');
        elseif((m>=1)&&(n>=1))
            enot2 = totg/(-n*F);
            ctemp2 = c1/n;
            pH2 = 0:1:14;
            tote2 = enot2-(ctemp2*totlog)-(ctemp2*m*pH2)-(ctemp2*tppog);
            EwpH(f3,:) = tote2;
            f3=f3+1;
            plot(pH2, tote2,'-m');
        end
    end
end
hold off

%Boundaries for number of equations = 5
if(noe<=5)
    figure('Name','Pourbaix Diagram','NumberTitle','off'); 
    hold on
    title(['Pourbaix Diagram at ',num2str(temp),' K, ',num2str(Conc),' M, ',num2str(ppog),' atm'])
    axis([dp0(1) dp0(2) dp0(3) dp0(4)])
    pHx = 0:1:14;
    E_H2 = -(0.0591*pHx);
    E_O2 = 1.23-(0.0591*pHx);
    plot(pHx,E_H2,'--.r',pHx,E_O2,'--.r');
    xlabel('Potential of Hydrogen (pH) ----->');
    ylabel('Potential (V) ----->');
    if(horil==1 && vertil==2 && slopl==2)
        minpH = min(pH);
        maxpH = max(pH);
        h1 = 0:0.0001:minpH;
        yh1 = Ea*ones(size(h1));
        plot(h1,yh1,'LineWidth',1.5);
        sl1 = EwpH(1,:);
        sl2 = EwpH(2,:);
        md1 = linspace(minpH,maxpH,length(sl1));
        md2 = linspace(maxpH,14,length(sl2));
        l1 = EwpH(1,pHx >= minpH & pHx <= maxpH);
        l2 = EwpH(1,pHx >= maxpH & pHx <= 14);
        l3 = EwpH(2,pHx >= minpH & pHx <= maxpH);
        l4 = EwpH(2,pHx >= maxpH & pHx <= 14);
        esl1 = linspace(l1(1),l1(length(l1)),15);
        esl2 = linspace(l2(1),l2(length(l2)),15);
        esl3 = linspace(l3(1),l3(length(l3)),15);
        esl4 = linspace(l4(1),l4(length(l4)),15);
        s1 = diff(sl1)/diff(md1);
        s2 = diff(sl2)/diff(md1);
        [x1,y1] = polyxpoly(pHx,sl1,pHx,sl2);
        if s1 > s2
            plot(md2,esl4,'LineWidth',1.5);
            plot(md1,esl1,'LineWidth',1.5);
        else 
            plot(md1,esl2,'LineWidth',1.5);
            plot(md2,esl3,'LineWidth',1.5);
        end
        v1 = dp0(4):-0.0001:Ea;
        xv1 = minpH*ones(size(v1));
        plot(xv1, v1,'LineWidth',1.5);
        v2 = dp0(4):-0.0001:y1;
        xv2 = x1*ones(size(v2));
        plot(xv2, v2,'LineWidth',1.5);        
    else
        fprintf('\nPlease Enter the Equation with 1 Horizontal, 2 Vertical and 2 Sloped Lines');
    end
else
    fprintf('\nStill Working on Equations > 5!');
end
