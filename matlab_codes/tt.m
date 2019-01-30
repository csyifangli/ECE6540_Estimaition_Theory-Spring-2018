% create data
    for ik = 1:4
        x(:,ik) = (-2:0.01:2)*pi/4*ik;
        A(:,ik) = sin(x(:,ik) + ik*pi/4);
    end
    % plot data
    fh = figure;
    ah = axes('Parent',fh);
    hold on
    for ik = 1:4
        plot(ah,x(:,ik),A(:,ik))
    end