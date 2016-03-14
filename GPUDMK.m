function Xkf= DMKalman(Vref,connectivity,transfer_matrix,observations,states,...
    excitable_threshold,exciting_threshold,deltaIndex,lambda)
%in this versio of DMK algotirhm I will put following variables into
%gpu:---they are transfer matrix,observation, first couple of state of
%states,Q,R, A, B
g=gpuDevice(1);
wait(g);
nodeNum = size(transfer_matrix,1);
%w = sqrt(nodeNum);
cellNum = size(transfer_matrix,2);
time=size(states,2); %the length of state sequence over time
%estimation=zeros(cellNum,time);
%Q = zeros(cellNum,cellNum);%process noise is 0

%R=zeros(nodeNum,nodeNum);%obsevation noise is zero
%W=0;
%V=zeros(nodeNum,1);%observation erro will be zero for now
%process model
%A = eye(cellNum); %process matirx
%B = eye(cellNum)'; % input control matrix
%U = zeros(cellNum,1); % no stimulus from outside, if it does, it will be in the form of state vector
%H = gpuArray(transfer_matrix);
%X = zeros(cellNum,time);% true state stroage sequence
%X = states;%assign the state sequence to the model
%X(:,1) = reshape(reshape([ones(w),ones(cellNum-w)],w,w),cellNum,1);%initialize the tissue

%P0 = zeros(cellNum,cellNum);% covriance matrix
%Z = observations;%the observation sequence
%Z(:,1)= H*X(:,1); %initialize observation
%Xkf = zeros(cellNum,time,'gpuArray');% state estimate initialize
%Xkf(:,1) = zeros(size(X(:,1)));
Xkf(:,1) = states(:,1);
Xkf(:,2) = states(:,2);
% err_P = zeros(time,cellNum);
% err_P(1,:) = (P0*ones(cellNum,1));%initialize the error overtime
%g.AvailableMemory
%I = eye(cellNum,'gpuArray');
%covA =1;
%randomA=sqrt(covA)*randn(2,2,time);
%covB =1;
%randomB=sqrt(covB)*randn(2,1,time);
%Gains = zeros(2,time);

for k=3:time
    Q = 100*diag(ones(1,cellNum,'gpuArray'));
    R = eye(nodeNum,'gpuArray');
    H = gpuArray(transfer_matrix);
    P0 = Q;
    I = eye(cellNum);
    A = diag(getDropRate(Vref,Xkf(:,k-1),deltaIndex));
    
    B = getControlMatrix(Vref,connectivity,Xkf(:,k-1),deltaIndex,excitable_threshold,exciting_threshold);
    A = gpuArray(A);
    
    B = gpuArray(B);
    gXkf=gpuArray(Xkf(:,k-1));
    Z=gpuArray(observations(:,k));
    
    %Z(:,k)= H*X(:,k)+V(k);%ovservation vector for one time instance
    %kalman filtering
    X_pre=A*gXkf+B*gXkf;% prediction of state
    P_pre = A*P0*A'+Q;% prediction of covriance
    Kg = P_pre*H'*inv(H*P_pre*H'+R);%kalman gain last term Sk
    %Kg = (H*P_pre*H'+R)*(H*P_pre*H'+R)'+lambda*lambda*...
    %        inv((Z(:,k)-H*X_pre)*(Z(:,k)-H*X_pre)')*(H*P_pre*H'+R)*H*P_pre;
    %regularised gain setting 1
%     Kg = (((H*P_pre*H'+R)*(H*P_pre*H'+R)'+(lambda*lambda).*...
%        ((Z(:,k)-H*X_pre)*(Z(:,k)-H*X_pre)'))\(H*P_pre*H'+R)*H*P_pre)';
    %reged gain setting 2
    %   Kg=zeros(cellNum,nodeNum);
    %gain setting 3
    
    %gain setting 4
    %Gains(:,k)=Kg;
    gXkf= X_pre+Kg*(Z-H*X_pre);%state innovation last term yk
    
    P0=(I-gather(Kg*H))*gather(P_pre);%covariance innovation
    %err_P(k,:) = (P0*ones(cellNum,1));%store the error covariance
    Xkf(:,k) = gather(gXkf);
    reset(g);
    wait(g);
    clc,disp(['generating estimation for ' num2str(k) ' msec']);
end

end
