function [err_P,measure,state_kf,Gains] = kf (A,B,U,transfer_matrix,state,Q,R,V)
nodeNum = size(transfer_matrix,1);
w = sqrt(nodeNum);
cellNum = size(state,1);
time=size(state,2); %the length of state sequence over time
%Q = zeros(cellNum,cellNum);%process noise is 0
%R=zeros(nodeNum,nodeNum);%obsevation noise is zero
%W=0;
%V=zeros(nodeNum,1);%observation erro will be zero for now
%process model
%A = eye(cellNum); %process matirx
%B = eye(cellNum)'; % input control matrix
%U = zeros(cellNum,1); % no stimulus from outside, if it does, it will be in the form of state vector
H = transfer_matrix;
%X = zeros(cellNum,time);% true state stroage sequence
X = state;%assign the state sequence to the model
%X(:,1) = reshape(reshape([ones(w),ones(cellNum-w)],w,w),cellNum,1);%initialize the tissue
P0 = [0 0;0 0];
%P0 = 0.25*eye(cellNum);% covriance matrix
Z = zeros(nodeNum,time);%the observation sequence
Z(:,1)= H*X(:,1); %initialize observation
Xkf = zeros(cellNum,time);% state estimate initialize
%Xkf(:,1) = zeros(size(X(:,1)));
Xkf(:,1) = state(:,1);
err_P = zeros(time,cellNum);
err_P(1,:) = (P0*ones(cellNum,1));%initialize the error overtime
I = eye(cellNum);
covA =1;
randomA=sqrt(covA)*randn(2,2,time);
covB =1;
randomB=sqrt(covB)*randn(2,1,time);
Gains = zeros(2,time);
for k=2:time
    A = [1,1;0,1];
    B = [0.5;1];
    Z(:,k)= H*X(:,k)+V(k);%ovservation vector for one time instance
    %kalman filtering
    X_pre=A*Xkf(:,k-1)+B*U;% prediction of state
    P_pre = A*P0*A'+Q;% prediction of covriance
    Kg = P_pre*H'*inv(H*P_pre*H'+R);%kalman gain
    Gains(:,k)=Kg;
    Xkf(:,k)= X_pre+Kg*(Z(k)-H*X_pre);%state innovation
    P0=(I-Kg*H)*P_pre;%covariance innovation
    err_P(k,:) = (P0*ones(cellNum,1));%store the error covariance
end
save('gains.mat','Gains');
state_kf=Xkf;
measure = Z;
end
    
    
    

