clc
prompt = 'Temperature: ';
Temp = input(prompt);
prompt = 'Standard Enthaply at 298 K: ';
Sen = input(prompt);
prompt = 'Standard Entropy at 298 K: ';
Set = input(prompt);
prompt = 'Molar Heat Capacity if given 1 and if not 0: ';
chc = input(prompt);
if (Temp==298)
        Sge = Sen - (Temp*Set);
        sprintf('%8.f',Sge)
elseif ((chc==0) && (Temp>298 || Temp<298))
        prompt = 'A Value: ';
        a = input(prompt);
        prompt = 'B Value: ';
        b = input(prompt);
        prompt = 'C Value: ';
        c = input(prompt);
        prompt = 'D Value: ';
        d = input(prompt);
        fun1 = @(t) a+b*t+c*power(t,2)+d*power(t,-2);
        q = integral(fun1,298,Temp);
        fun2 = @(t) a*power(t,-1)+b+c*power(t,1)+d*power(t,-3);
        p = integral(fun2,298,Temp);
        Sge1 = (Sen + q);
        Sge2 = -(Temp*(Set+p));
        Sge = Sge1+Sge2;
        sprintf('%8.f',Sge)
elseif ((chc==1) && (Temp>298 || Temp<298))
        prompt = 'Molar Heat Capacity at 298 K: ';
        shc = input(prompt);
        fun1 = @(t) shc*power(t,0);
        q = integral(fun1,298,Temp);
        fun2 = @(t) shc*power(t,-1);
        p = integral(fun2,298,Temp);
        Sge1 = (Sen + q);
        Sge2 = -(Temp*(Set+p));
        Sge = Sge1+Sge2;
        sprintf('%8.f',Sge)        
end
