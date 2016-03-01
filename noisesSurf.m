Q1 = 0:1:100;
Q2 = 0:0.005:0.5; %speed cov zoomed  times
R=100;

pos_err=zeros(size(Q1,2),size(Q2,2));
speed_err=zeros(size(Q1,2),size(Q2,2));
pos_cov=zeros(size(Q1,2),size(Q2,2));
speed_cov=zeros(size(Q1,2),size(Q2,2));
for i=1:size(Q1,2)
    for j=1:size(Q2,2)
        %[pos_err(i,j),speed_err(i,j),pos_cov(i,j),speed_cov(i,j)]...
        [a,b,c,d]=freefall([Q1(i) 0;0 Q2(j)],R);
        pos_err(i,j)=a;%result(1,1);
        speed_err(i,j)=b;%result(1,2);
        pos_cov(i,j)=c;%result(1,3);
        speed_cov(i,j)=d;% result(1,4);
        
    end
end

% figure
% subplot(2,2,1)
% surf(pos_err)
% title('pos_err')
% subplot(2,2,2)
% surf(speed_err)
% title('speed_err')
% subplot(2,2,3)
% surf(pos_cov)
% title('pos_cov')
% subplot(2,2,4)
% surf(speed_cov)
% title('speed_cov')
% 




