
function [f0,g0,h0]= CH1
%Q1
    function f
        f(x) = sin(x-1);
    end
    function g
        g(x)=(x-1).^(2/3);
    end
    function H 
        if x<=0
            H(x) = x.^3;
        else
            H(x) = -x.^3;
        end
    end
    x = -10:1:10;
    %x
    arrf = zeros(1,length(x));
    arrh = zeros(1,length(x));
    for i=1:length(x)
        arrf(i) = sin(x(i)-1);
        if x(i)<=0
            arrh(i) = 2*x(i).^3;
        else
            arrh(i) = -x(i).^3;
        end
    end
    %plot(x,arrf)
    %plot(x,(x-1).^(2/3))
    %plot(x,arrh)
    %df = diff();
    %dg = diff(g);
    %dH = diff(H);
    %f0 = subs(df,x,0);
    %g0 = subs(dg,x,0);
    %h0 = subs(dH,x,0);
    
%Q2
    h = [0.1,0.05,0.01,0.005,0.001];
    f0=sin(-1);
    arr_err=[0,0,0,0,0];
    for i=1:length(h)
        %fhr = sin(h(i)-1); %f(x_0+h)
        %fhl = sin(-h(i)-1);%f(x_0-h)
        %fh = -h(i).^3;
        fhr = -h(i).^3;
        fhl = 2*h(i).^3;
        %df_apx = (fh-0)/h(i);
        df_apx = (fhr-fhl)/(2*h(i));
        arr_err(i) = abs(-df_apx);
    end
    arr_err
    plot(log(h),log(arr_err))
    %dg_apx = (g(h)-g(0))/h
    %dH_appr = (H(h)-H(0))/h
    
       
    
end
