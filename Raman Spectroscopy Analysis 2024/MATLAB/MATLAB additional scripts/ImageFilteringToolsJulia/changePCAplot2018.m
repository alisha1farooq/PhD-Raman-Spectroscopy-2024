function changePCAplot2018

global PCARun

a=gcf;PCAnumber=a.Children(2).Value;
disp(num2str(PCAnumber))
figure(gcf);
subplot(2,1,1)
cla
imagesc(PCARun.PCAScores(:,:,PCAnumber))
title(['PC Score: ',num2str(PCAnumber)])
axis image   
colorbar
subplot(2,1,2)
cla
plot(PCARun.xj,PCARun.PCALoadings(PCAnumber,:),'tag','PCloading')
hold on
plot(PCARun.xj,zeros(size(PCARun.xj,2)),'k')
hold off
title('PCloading')
end
