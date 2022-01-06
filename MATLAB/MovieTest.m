h = figure;
% The variable loops is the number of frames
loops = 20;
M(loops) = struct('cdata',[],'colormap',[]);
h.Visible = 'off';

for i = 1:20
    scale = 1 - 0.05*i;
    [X, Y] = meshgrid(1:0.5:10,1:20);
    Z = scale * (sin(X) + cos(Y));
    % Use surf for movie, since movie command works with colormap objects,
    % which are not associated with plot or scatter command.
    surf(X,Y,Z)
    % Make sure to set axes limits AFTER creating the surface plot
    % Axes limits must be standardized for the movie (otherwise the scales
    % in each frame will change.
    xlim([0,10])
    ylim([0,20])
    zlim([-2,2])
    drawnow
    M(i) = getframe;
end

h.Visible = 'on';
% Second argument = number of times to loop the movie
% Third argument = fps, so larger fps means faster animation
% Feel free to play around with these parameters to see how this movie
% command works
movie(M, 1, 10);