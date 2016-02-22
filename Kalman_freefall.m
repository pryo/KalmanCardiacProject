tic
N = 1000;
time = 1:1:N;
X =zeros(2,N);
seudoA = [0 1.5;1.5 0];
A = [1,1;0,1];%process model for freefall
X(:,1) = [95;1];%initial state height 95m, speed 1m/s
B=[0.5;1];
U=-1;
Q=[1,0;0,1]; %process noise covariance;
R = 1;%observation covariance
W = sqrt(Q)*randn(2,N);%process noise
V=sqrt(R)*randn(1,N);%observation noise
for k=2:N
    X(:,k)= A*X(:,k-1)+B*U+W(k);
end
H = [1,0];
[err_P,Z,Xkf] = kf(A,B,U,H,X,Q,R,V);

for k =1:N
    measure_err(k) = Z(k)-X(1,k);
    kalman_err_x(k)= Xkf(1,k)-X(1,k);
    kalman_err_v(k) = Xkf(2,k)-X(2,k);
end
toc
figure('name','process noise')
subplot(2,1,1)
plot(W(1,:));
xlabel('time/s');
ylabel('position noise');
subplot(2,1,2)
plot(W(2,:));
xlabel('time/s');
ylabel('speed noise');

figure('name','measurement noise')
plot(V);
xlabel('time/s');
ylabel('measure noise');

figure('name','kalman error')
subplot(2,1,1)
hold on,box on;
plot(measure_err,'-r.');%measurement error on position
plot(kalman_err_x,'-g.');%estimated position error
legend('measured position error','estimated position error')
xlabel('time/s');
ylabel('error/m');

subplot(2,1,2)
plot(kalman_err_v);
xlabel('time/s');
ylabel('speed error/m/s');

figure('name','state covariance')
subplot(2,1,1)
plot(err_P(:,1));
xlabel('time/s');
ylabel('position error covariance');

subplot(2,1,2)
plot(err_P(:,2));
xlabel('time/s');
ylabel('speed error covariance');



% figure('name','position error')
% plot(kalman_err_x);
% figure('name','speed error')
% plot(kalman_err_v);
% figure('name','position')
% plot(time,Xkf(1,:),'r',time,X(1,:),'g');
% figure('name','speed')
% plot(time,X(2,:),'g');
% hold on;
% plot(time,Xkf(2,:),'r');