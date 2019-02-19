clear all;
close all;
lJV = @(r, sigma,epsilon) (4.*epsilon.*((sigma./r).^12-(sigma./r).^6));
lJF = @(r, sigma,epsilon) (-4.*epsilon.*(-12.*(sigma.^12./r.^13)+6.*(sigma.^6./r.^7)));
rad = 0:1:1000;
plotPotential = [];
potential = [];
scalar = length(rad)/3;
for r_ab = rad;
    for r_cb = rad;
        potential(r_ab+1,r_cb+1) = lJV(r_ab/scalar,1,1)+lJV(r_cb/scalar,1,1);
        if potential(r_ab+1,r_cb+1)>0
            plotPotential(r_ab+1,r_cb+1)=NaN;
        else
            plotPotential(r_ab+1,r_cb+1)=potential(r_ab+1,r_cb+1);
        end
    end
end

f=surf(rad/scalar,rad/scalar,plotPotential);
set(f,'LineStyle','none');
set(gcf, 'Position',  [0, 0, 900, 900]);
view(0,90);
drawnow;
hold on;

m_a = 1;
m_b = 1;
m_c = 1;

x_a = -1.3;
x_b = 0;
x_c = 2.5;
v_a = 0;
v_b = 0;
v_c = 0;

r_abList = [];
r_cbList = [];
rv_abList = [];
rv_abMax = 0;
rv_cbList = [];
rv_cbMax = 0;
x_aList = [];
x_bList = [];
x_cList = [];
timeStep = 0.001;
timeLimit = 9;
time = 0:timeStep:timeLimit;

for t = time
    F_a= -lJF(abs(x_b-x_a),1,1)-lJF(abs(x_c-x_a),1,1);
    F_b= lJF(abs(x_a-x_b),1,1)-lJF(abs(x_c-x_b),1,1);
    F_c= lJF(abs(x_b-x_c),1,1)+lJF(abs(x_c-x_a),1,1);
    
    x_a = x_a + v_a*timeStep + F_a/m_a*timeStep^2/2;
    x_b = x_b + v_b*timeStep + F_b/m_b*timeStep^2/2;
    x_c = x_c + v_c*timeStep + F_c/m_c*timeStep^2/2;
    
    v_a = v_a+(F_a+(-lJF(abs(x_b-x_a),1,1)-lJF(abs(x_c-x_a),1,1)))/m_a*timeStep/2;
    v_b = v_b+(F_b+(lJF(abs(x_a-x_b),1,1)-lJF(abs(x_c-x_b),1,1)))/m_b*timeStep/2;
    v_c = v_c+(F_c+(lJF(abs(x_b-x_c),1,1)+lJF(abs(x_c-x_a),1,1)))/m_c*timeStep/2;
    
    r_abList = [r_abList abs(x_b-x_a)];
    r_cbList = [r_cbList abs(x_c-x_b)];
    
    rv_abList = [rv_abList abs(v_b-v_a)];
    rv_cbList = [rv_cbList abs(v_c-v_b)];
    
    if rv_abList(end) > rv_abMax
        rv_abMax = rv_abList(end);
    end
    
    if rv_cbList(end) > rv_cbMax
         rv_cbMax = rv_cbList(end);
    end
    
    x_aList = [x_aList x_a];
    x_bList = [x_bList x_b];
    x_cList = [x_cList x_c];
end

for i = 1:1:length(r_abList)-1
    line([r_abList(i) r_abList(i+1)],[r_cbList(i) r_cbList(i+1)],[(lJV(r_abList(i),1,1)+lJV(r_cbList(i),1,1)) (lJV(r_abList(i+1),1,1)+lJV(r_cbList(i+1),1,1))],'color',[rv_cbList(i)/rv_cbMax rv_abList(i)/rv_abMax 0])
end

%plot(r_abList,r_cbList);
figure();
hold on;
c = zeros(3,length(x_aList));
c = c+1;
plot(time,x_aList);
plot(time,x_bList);
plot(time,x_cList);










