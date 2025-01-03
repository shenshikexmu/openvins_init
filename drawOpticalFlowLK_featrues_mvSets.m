function  drawOpticalFlowLK_featrues_mvSets(imgpyr,frame1,frame2,pts1,pts2,mvSets)


img1=imgpyr{frame1}{1};
img2=imgpyr{frame2}{1};


colors = jet(6);

figure;
imshow([img1,img2;img2,img2]);
hold on;

s=0;
for m=1:size(mvSets,1)  

    k=mvSets(m);
 
    s=s+1;
    hold on
    plot(pts1(k,1),pts1(k,2),'r+');
    
    if mod(s,2)==0
        hold on;
        plot(pts2(k,1),pts2(k,2)+size(img1,1),'g+');
        hold on;
        plot([pts1(k,1),pts2(k,1)],[pts1(k,2),pts2(k,2)+size(img1,1)],'-', 'Color', colors(mod(s,6)+1,:));
    else
        hold on;
        plot(pts2(k,1)+size(img1,2),pts2(k,2),'g+');
        hold on;
        plot([pts1(k,1),pts2(k,1)+size(img1,2)],[pts1(k,2),pts2(k,2)],'-', 'Color', colors(mod(s,6)+1,:));
    end

    hold on;
    %quiver(pts{1}(j,1)+size(img1,2), pts{1}(j,2)+size(img1,1), pts{n}(k,1)-pts{1}(j,1), pts{n}(k,2)-pts{1}(j,2), 'r', 'LineWidth', 2);

    plot([pts1(k,1)+size(img1,2),pts2(k,1)+size(img1,2)],[pts1(k,2)+size(img1,1),pts2(k,2)+size(img1,1)], 'r', 'LineWidth', 2);
    
end

title(['Optical Flow (Lucas-Kanade) image',num2str(frame1),' to image',num2str(frame2)]);

hold off;


end
    