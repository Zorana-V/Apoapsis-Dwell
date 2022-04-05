function mercatorDisplay(lonE,lat)
earth = imread('earth.jpg');
% Open a new figure, then run the command "clf"
% Example:
% figure(<next-figure-number-here>);
% clf
% Then run the code below.
image('CData',earth,'XData',[-180 180],'YData',[90 -90])
hold on
plot(lonE*180/pi,lat*180/pi,'w*');
end