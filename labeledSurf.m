figure
subplot(2,2,1)
surf(pos_err)
xlabel('speed error cov')
ylabel('position error cov')
title('pos_err')
subplot(2,2,2)
surf(speed_err)
xlabel('speed error cov')
ylabel('position error cov')
title('speed_err')
subplot(2,2,3)
surf(pos_cov)
xlabel('speed error cov')
ylabel('position error cov')
title('pos_cov')
subplot(2,2,4)
surf(speed_cov)
xlabel('speed error cov')
ylabel('position error cov')
title('speed_cov')