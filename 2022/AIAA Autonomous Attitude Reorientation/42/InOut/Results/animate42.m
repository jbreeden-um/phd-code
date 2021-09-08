% Use the CamSnap feature to record the images necessary for this file to run 
t0 = 0;
t1 = 2100;
frames_per_simtime = 1;
video_duration = 100;
nframes = (t1-t0)*frames_per_simtime;
frame_rate = nframes/video_duration;

record = 1;
data = load('results42.mat');

moviename = 'MyMovie1';
vidObj = VideoWriter(moviename);
vidObj.FrameRate = frame_rate;
if record, open(vidObj); end

fig = figure(1); clf;
set(fig, 'Position', [1000 60 1600 800]);

im1 = imread(sprintf('42 without QP/CamFrame%05d.ppm',1));
im2 = imread(sprintf('42 with QP/CamFrame%05d.ppm',1));
im = imfuse(im1, im2, 'montage');
im(:) = 0;
imshow(im);
text(40, 150, ['"Autonomous Spacecraft Attitude Reorientation' newline ...
    '                Using Robust Sampled-Data' newline ...
    '                  Control Barrier Functions"'],...
    'Color', 'w', 'Fontsize', 50, 'FontWeight', 'Bold');
text(290, 480, 'Comparison in 42 Spacecraft Simulator', 'Color', 'w', 'Fontsize', 44);
text(297, 488, '_______________________________', 'color', 'w', 'fontsize', 44);
text(460, 650, 'Joseph Breeden and Dimitra Panagou', 'Color', 'w', 'Fontsize', 30);

duration = 5;
fade = 1;
if record
    for i=1:duration*frame_rate
        this_frame = getframe(fig);
        this_frame.cdata = this_frame.cdata(30:829,87:1686,:);
        writeVideo(vidObj,this_frame);
    end
    fade_frames = round(fade*frame_rate);
    for i=1:fade_frames
        this_frame = getframe(fig);
        new_im = this_frame.cdata(30:829,87:1686,:);
        m = max(new_im(:));
        new_im = new_im - (m*i/fade_frames);
        imshow(new_im);
        this_frame = getframe(fig);
        this_frame.cdata = this_frame.cdata(30:829,87:1686,:);
        writeVideo(vidObj,this_frame);
    end
end

imshow(im);
% The yellow vector is the current Sun vector. The pale blue vector is the
% current Moon vector. The instrument must avoid pointing at the Sun by 25^o, and
% the star tracker must avoid pointing at the Sun by 45^o and the Moon by
% 30^o. The following is a comparison with and without using ZohCBFs.
text(140,150,'The yellow vector is the current Sun vector.','Color',[1;1;0],'FontSize',40);
text(140,200,'The pale blue vector is the current Moon vector.','Color',[50;50;100]/150,'FontSize',40);
text(140,270,'The instrument must avoid pointing at the Sun by 25^o.','color','w','fontsize',40); 
text(140, 335, 'The star tracker must avoid pointing at the Sun','color','w','fontsize',40); 
text(190, 390, 'by 45^o and the Moon by 30^o.','color','w','fontsize',40); 
text(140, 480, 'The red vectors are an inertial frame.','color',[0.8;0;0],'fontsize',40);
text(140, 580, 'The following is a comparison with and without','color','w','fontsize',40);
text(190, 635, 'using ZohCBFs for safety.','color','w','fontsize',40); 

duration = 6;
fade = 1;
if record
    for i=1:duration*frame_rate
        this_frame = getframe(fig);
        this_frame.cdata = this_frame.cdata(30:829,87:1686,:);
        writeVideo(vidObj,this_frame);
    end
    fade_frames = round(fade*frame_rate);
    for i=1:fade_frames
        this_frame = getframe(fig);
        new_im = this_frame.cdata(30:829,87:1686,:);
        m = max(new_im(:));
        new_im = new_im - (m*i/fade_frames);
        imshow(new_im);
        this_frame = getframe(fig);
        this_frame.cdata = this_frame.cdata(30:829,87:1686,:);
        writeVideo(vidObj,this_frame);
    end
end

% The following is to help with alignment.
% hold on;
% for i=0:100:1600
%    plot([i i], [0 800], 'w');
% end

% nframes = 10;
for frame=6:nframes
    im1 = imread(sprintf('42 without QP/CamFrame%05d.ppm',frame));
    im2 = imread(sprintf('42 with QP/CamFrame%05d.ppm',frame));
    im = imfuse(im1, im2, 'montage');
    imshow(im);
    text(280, 40, 'Shortest Path', 'Color', 'w', 'FontSize', 30);
    text(1100, 40, 'ZohCBF-QP', 'Color', 'w', 'FontSize', 30);
    
    % Nominal
    h = data.h_nom(frame,1); if h <= 0, color = [0 1 0.8]; else, color = 'r'; end
    d = acosd(h + cosd(25));
    text(5, 710, sprintf('Instrument to Sun: %.2f^o', d), 'Color', color, 'FontSize', 14)
    
    h = data.h_nom(frame,2); if h <= 0, color = [0 1 0.8]; else, color = 'r'; end
    d = acosd(h + cosd(45));
    text(5, 735, sprintf('Star Tracker to Sun: %.2f^o', d), 'Color', color, 'FontSize', 14)
    
    h = data.h_nom(frame,3); if h <= 0, color = [0 1 0.8]; else, color = 'r'; end
    d = acosd(h + cosd(30));
    text(5, 760, sprintf('Star Tracker to Moon: %.2f^o', d), 'Color', color, 'FontSize', 14)
    
    h = data.Safety_nom(frame,7); if h <= 0, color = [0 1 0.8]; else, color = 'r'; end
    text(5, 785, sprintf('Kinetic Energy: %.4f mJ', (h+5.1e-05)*1e3), 'Color', color, 'FontSize', 14)
    
    % ZohCBF
    h = data.h_zoh(frame,1); if h <= 0, color = [0 1 0.8]; else, color = 'r'; end
    d = acosd(h + cosd(25));
    text(805, 710, sprintf('Instrument to Sun: %.2f^o', d), 'Color', color, 'FontSize', 14)
    
    h = data.h_zoh(frame,2); if h <= 0, color = [0 1 0.8]; else, color = 'r'; end
    d = acosd(h + cosd(45));
    text(805, 735, sprintf('Star Tracker to Sun: %.2f^o', d), 'Color', color, 'FontSize', 14)
    
    h = data.h_zoh(frame,3); if h <= 0, color = [0 1 0.8]; else, color = 'r'; end
    d = acosd(h + cosd(30));
    text(805, 760, sprintf('Star Tracker to Moon: %.2f^o', d), 'Color', color, 'FontSize', 14)
    
    h = data.Safety_zoh(frame,7); if h <= 0, color = [0 1 0.8]; else, color = 'r'; end
    text(805, 785, sprintf('Kinetic Energy: %.4f mJ', (h+5.1e-05)*1e3), 'Color', color, 'FontSize', 14)
    
    drawnow;
    
	if record
        this_frame = getframe(fig);
        this_frame.cdata = this_frame.cdata(30:829,87:1686,:);
        writeVideo(vidObj,this_frame);
    end
    
    waitbar(frame/nframes);
end

imshow(im*0);
text(150, 400, 'Thank you to NASA for the simulation and graphics engine.', 'color', 'w', 'fontsize', 36);

duration = 2;
fade = 0.5;
if record
    for i=1:duration*frame_rate
        this_frame = getframe(fig);
        this_frame.cdata = this_frame.cdata(30:829,87:1686,:);
        writeVideo(vidObj,this_frame);
    end
    fade_frames = round(fade*frame_rate);
    for i=1:fade_frames
        this_frame = getframe(fig);
        new_im = this_frame.cdata(30:829,87:1686,:);
        m = max(new_im(:));
        new_im = new_im - (m*i/fade_frames);
        imshow(new_im);
        this_frame = getframe(fig);
        this_frame.cdata = this_frame.cdata(30:829,87:1686,:);
        writeVideo(vidObj,this_frame);
    end
end

if record, close(vidObj); end