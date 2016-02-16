N = 1000;
time = 1:1:N;
X =zeros(2,N);
A = [1,1;0,1];%process model for freefall
X(:,1) = [95;1];%initial state height 95m, speed 1m/s
B=[0.5;1];
U=-1;
Q=[2,2;2,2]; %process noise covariance;
R = 1;
W = sqrt(Q)*randn(2,N);
V=sqrt(R)*randn(1,N);%observation noise
for k=2:N
    X(:,k)= A*X(:,k-1)+B*U+W(k);
end
H = [1,0];
Xkf = kf(A,B,U,H,X,Q,R,V);

for k =1:N
    kalman_err_x(k)= Xkf(1,k)-X(1,k);
    kalman_err_v(k) = Xkf(2,k)-X(2,k);
end
% figure
% plot(W);
% xlabel('time/s');
% ylabel('process noise');
% figure
% plot(V);
% xlabel('time/s');
% ylabel('measure noise');
% figure
% hold on,box on;
figure('name','position error')
plot(kalman_err_x);
figure('name','speed error')
plot(kalman_err_v);
figure('name','position')
plot(time,Xkf(1,:),time,X(1,:),'g');
figure('name','speed')
plot(time,X(2,:),'g');
hold on;
plot(time,Xkf(2,:),'r');