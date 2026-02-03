
% load the datafile scatterPlotFits.mat & pennyColor.mat 
j = 3 ; % 1=apt, 2=zeb, 3 =cfx, 4 =cfc
figure();
for i=1:3
    transp = [1, 1,1];
    markers = ['o','^','s'];

    for ii=1:4
     
%         ss(ii,i) = scatter3(c(ii,1,i,j),b(ii,1,i,j),a(ii,1,i,j),200,cly(ii,:,i),'filled','Marker', markers(i));
%         set(ss(ii,i), 'MarkerFaceAlpha', transp(i))
%         hold on; %if not working for the colors, use ncly 
%         
        ss(ii,i) = scatter3(cx(ii,1,i,j),bx(ii,1,i,j),ax(ii,1,i,j),200,clx(ii,:),'filled','Marker', markers(i));
        set(ss(ii,i), 'MarkerFaceAlpha', transp(i))
        hold on;
        formatFigure('ccg width', 'Latency', '');
         xlabel('ccg width','FontSize',35,'Color','k','FontWeight','bold'); 
         ylabel('Latency','FontSize',35,'Color','k','FontWeight','bold');
         zlabel('Peak','FontSize',35,'Color','k','FontWeight','bold')
       
    end;
end;

% for i=1:4
%     transp = [1, 1,1];
%     markers = ['o','^','s'];
% 
%     for ii=1:3
%         
%         ss(i,ii) = scatter3(cx(i,1,ii,j),bx(i,1,ii,j),ax(i,1,ii,j),100,clx(i,:),'filled','Marker', markers(ii));
%         set(ss(i,ii), 'MarkerFaceAlpha', transp(ii))
%         hold on;
%          xlabel('ccg width','FontSize',16); ylabel('Latency','FontSize',16);zlabel('Peak','FontSize',16)
%        
%     end;
% end;