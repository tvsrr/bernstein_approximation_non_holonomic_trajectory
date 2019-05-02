%Bernstein's path creation 
%Author- Raghu Ram Theerthala 
%date- 28/04/19
%Matlab verison - R2018b

% declare the initial state of the vehicle 
x0 = 0; y0 = 0; th0 = 0; t0 = 0; xdot0=0;
xf = 10; yf = 10; thf = 0; tf = 10; xdotf=2;
kdot0 = 0;ytc = 0; ydotc=1; 
kdotf = 0.5;xtc = 5; xdotc=1; tc=5;

%ydot(t) =  xdot(t)tan(theta(t)) is the non holonomic condition 
%x(t) and tan(theta(t)) are non integrabale. Thus we need an approximation
%We choose bernsteins approximation here. 
order = 5; 

%declare weights 
x_weights = sym('Wx',[order+1 1]);
theta_weights = sym('Wt',[order+1 1]);

%declare bernstein function 
syms mu(t) B(l)
mu(t) = (t-t0)/(tf-t0);
B(l) = nchoosek(order,l)*((1-mu)^(order-l))*mu^(l);
Bdot(l) = diff(B(l),t);

%getting xdot(t)
syms xdot(t) x(t) tan_approx(t) tandot_approx(t)
k=0;fi=0;
for i = 0:order
k = k + x_weights(i+1)*(B(i));
fi=fi+x_weights(i+1)*Bdot(i);
end
x(t)=k;xdot(t)=fi;

%getting tan(theta(t)) approximation
ta = 0;tadot = 0;
for i =0:order
ta=ta+theta_weights(i+1)*B(i);
tadot = tadot + theta_weights(i+1)*Bdot(i);
end
tan_approx(t)=ta;
tandot_approx(t)=tadot;


%We need to add constraints to get weights
%say for order 5 
%You can add constraints here.
%This code is solved for 5th order

%% Finding Wx weights 
%determine Wx1 to Wx6 
x_weights(1)=x0;
x_weights(6)=xf;

R1(t)=xtc-x_weights(1)*B(0)-x_weights(6)*B(5);
R2(t)=xdot0-x_weights(1)*Bdot(0)-x_weights(6)*Bdot(5);
R3(t)=xdotf-x_weights(1)*Bdot(0)-x_weights(6)*Bdot(5);
R4(t)=xdotc-x_weights(1)*Bdot(0)-x_weights(6)*Bdot(5);


Bx1(t)=[B(1) B(2) B(3) B(4)];
Bx2(t)=[Bdot(1) Bdot(2) Bdot(3) Bdot(4)];
Bx3(t)=Bx2(t);
Bx4(t)=Bx2(t);

Ax=[R1(tc);R2(t0);R3(tf);R4(tc)];
Bx=[Bx1(tc);Bx2(t0);Bx3(tf);Bx4(tc)];

Wx= pinv(Bx)*Ax;

x_weights(2:5)=Wx;
%% Finding Wk weights



fi=0;
for i = 0:order
k = k + x_weights(i+1)*(B(i));
fi=fi+x_weights(i+1)*Bdot(i);
end
xdot(t)=fi;

fun =(xdot*tan_approx);
% y(t) = y0 +int(fun,t,t0,t);
syms Fu(a1,a2,a3,a4,a5,a6) 

Fu(theta_weights(1),theta_weights(2),theta_weights(3),theta_weights(4),theta_weights(5),theta_weights(6)) = fun(t);
f1(t) = Fu(1,0,0,0,0,0);
F1(t) = int(f1(t),t,t0,t);

f2(t) = Fu(0,1,0,0,0,0);
F2(t) = int(f2(t),t,t0,t);

f3(t) = Fu(0,0,1,0,0,0);
F3(t) = int(f3(t),t,t0,t);

f4(t) = Fu(0,0,0,1,0,0);
F4(t) = int(f4(t),t,t0,t);

f5(t) = Fu(0,0,0,0,1,0);
F5(t) = int(f5(t),t,t0,t);

f6(t) = Fu(0,0,0,0,0,1);
F6(t) = int(f6(t),t,t0,t);

theta_weights(1)= tan(th0);
theta_weights(6)= tan(thf);

ta = 0;
for i =0:order
ta=ta+theta_weights(i+1)*B(i);
end
tan_approx(t)=ta;
tandot_approx(t)=diff(tan_approx(t),t);
Rk1 = [kdot0-theta_weights(1)*Bdot(0)-theta_weights(6)*Bdot(5)];
Rk2 = [kdotf-theta_weights(1)*Bdot(0)-theta_weights(6)*Bdot(5)];
Rk3 = ytc-theta_weights(1)*F1(tc)-theta_weights(6)*F6(tc);
Rk4 = yf-theta_weights(1)*F1(tf)-theta_weights(6)*F6(tf);


Ak = [Rk1;Rk2;Rk3;Rk4];

Bk1(t) = [Bdot(1) Bdot(2) Bdot(3) Bdot(4)];
Bk2(t) = Bk1(t);
Bk3(t) = [F2 F3 F4 F5];
Bk4(t) = Bk3(t);

Bk = [Bk1(t0);Bk2(tf);Bk3(tc);Bk4(tf)];

Wk = pinv(Bk)*Ak;

theta_weights(2:5)=Wk;

%% Plot section 
syms l(t) xdot_f(t) tant(t)
temp_x =0;temp_dotx=0;
temp_y = 0; i=1;
for k = 0:order
    temp_x = temp_x+x_weights(k+1)*B(k);
    temp_dotx = temp_dotx+x_weights(k+1)*Bdot(k);
    temp_y = temp_y+theta_weights(k+1)*B(k);
end
l(t)=temp_x; xdot_f(t)= temp_dotx; tant(t) = temp_y;
for ti=0:0.1:tf
x_f(i) = l(ti); 
y_f(i) = y0 + int(xdot_f(t)*tant(t),t,t0,ti);
i=i+1;
end
plot(x_f,y_f);
hold on
plot(x0,y0,'m.','markersize',30);
plot(xtc,ytc,'m.','markersize',30);
plot(xf,yf,'m.','markersize',30);









