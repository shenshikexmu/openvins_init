function  drawOpticalFlowLK_featrues(imgpyr,features,table_img_timestamp,cam_id1,cam_id2,frame1,frame2)


img1=imgpyr{frame1}{1};
img2=imgpyr{frame2}{1};

timestamp1=table_img_timestamp(frame1,2);
timestamp2=table_img_timestamp(frame2,2);

allKeys = keys(features);
k=0;
pts1=[];
pts2=[];
for i = 1:length(allKeys)

    key = allKeys{i}; % 获取当前键

    feat = features(key);

    frame1_timestamps=0;

    if size(feat.timestamps{cam_id1},1)>1

        for j=1:size(feat.timestamps{cam_id1},1)

            if timestamp1==feat.timestamps{cam_id1}(j)
                
                frame1_timestamps=j;

            end

        end

    end

    frame2_timestamps=0;

    if size(feat.timestamps{cam_id2},1)>1

        for j=1:size(feat.timestamps{cam_id2},1)

            if timestamp2==feat.timestamps{cam_id2}(j)
                
                frame2_timestamps=j;

            end

        end

    end

    if frame1_timestamps~=0 && frame2_timestamps~=0

        k=k+1;

        pts1=[pts1;feat.uvs{cam_id1}(frame1_timestamps,:)];

        pts2=[pts2;feat.uvs{cam_id2}(frame2_timestamps,:)];

    end

end


colors = jet(6);

figure;
imshow([img1,img2;img2,img2]);
hold on;

s=0;
for k=1:size(pts1,1)   
 
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

a=10;



end
    