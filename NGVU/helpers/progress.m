function status = progress(t, y, flag, varargin)
% Minimal ODE progress function (horizontal bar only)
% - One slim horizontal bar with percent, elapsed, and ETA
% - No nested helper functions
% - Throttled UI updates (~0.5 s)
%
% Use with: options = odeset('OutputFcn', @odeprog, ...);

persistent H   % holds UI handles/state
status = 0;    % never stop solver

% -------------- INIT / DONE --------------
if nargin >= 3 && ~isempty(flag)
    switch flag
        case 'init'
            % Final time for this ODE call (best available guess)
            if isempty(t), tfin = 1; else, tfin = t(end); end
            if tfin <= 0, tfin = 1; end

            % Create UI
            H = struct();
            H.tfin = tfin;
            H.tstart = clock;
            H.lastUpdate = 0;

            H.fig = figure(95); clf(H.fig);
            set(H.fig, 'Name','ODE Progress', 'NumberTitle','off', ...
                'Color','w', 'MenuBar','none', 'ToolBar','none', ...
                'Position',[200 200 600 120]);

            H.ax = axes('Parent',H.fig, 'Position',[0.07 0.35 0.86 0.45]);
            axis(H.ax, [0 1 0 1]);
            set(H.ax, 'XTick',[], 'YTick',[], 'Box','on');
            title(H.ax, '0%');

            % Background and foreground patches
            H.bg = patch('Parent',H.ax, 'XData',[0 1 1 0], 'YData',[0 0 1 1], ...
                         'FaceColor',[0.92 0.92 0.92], 'EdgeColor',[0.8 0.8 0.8]);
            H.fg = patch('Parent',H.ax, 'XData',[0 0 0 0], 'YData',[0 0 1 1], ...
                         'FaceColor',[0.30 0.60 1.00], 'EdgeColor','none');

            H.txt = uicontrol('Style','text','Parent',H.fig, ...
                'Units','normalized','Position',[0.07 0.1 0.86 0.18], ...
                'BackgroundColor','w','HorizontalAlignment','center', ...
                'String','Elapsed: 0s    |    ETA: --');

            drawnow;

        case 'done'
            % Final refresh and close
            if ~isempty(H) && isfield(H,'fig') && ishandle(H.fig)
                % Force 100% bar
                if isfield(H,'fg') && ishandle(H.fg)
                    set(H.fg,'XData',[0 1 1 0]);
                    set(get(H.ax,'Title'),'String','100%');
                end
                % Update elapsed/ETA text
                telapsed = etime(clock, H.tstart);
                % ETA = 0 at completion
                set(H.txt, 'String', sprintf('Elapsed: %s    |    ETA: 0s', local_fmt_time(telapsed)));
                drawnow;
                close(H.fig);
            end
            H = [];

        otherwise
            % ignore other flags
    end
    return;
end

% -------------- REGULAR UPDATES --------------
if isempty(H) || ~isfield(H,'fig') || ~ishandle(H.fig)
    return; % UI not initialized
end

% Throttle UI updates to ~0.5 s
if etime(clock, H.lastUpdate) < 0.5
    return;
end

% Current progress fraction
tnow = t(end);
pct = tnow / H.tfin;
if ~isfinite(pct), pct = 0; end
pct = min(max(pct,0),1);

% Update bar
if ishandle(H.fg)
    set(H.fg, 'XData', [0 pct pct 0]);
end
if ishandle(H.ax)
    set(get(H.ax,'Title'),'String', sprintf('%.1f%%', 100*pct));
end

% Elapsed + ETA
telapsed = etime(clock, H.tstart);
if pct > 0
    tremain = telapsed * (1/pct - 1);
else
    tremain = inf;
end
if ishandle(H.txt)
    set(H.txt,'String', sprintf('Elapsed: %s    |    ETA: %s', ...
        local_fmt_time(telapsed), local_fmt_time(tremain)));
end

H.lastUpdate = clock;  % persist last update time
drawnow limitrate;

end % function odeprog

% --------- Inline time formatter ---------
function s = local_fmt_time(seconds)
seconds = max(0, round(seconds));
dd = floor(seconds/86400); seconds = seconds - dd*86400;
hh = floor(seconds/3600);  seconds = seconds - hh*3600;
mm = floor(seconds/60);    ss = seconds - mm*60;
if dd > 0
    s = sprintf('%dd %02dh %02dm %02ds', dd, hh, mm, ss);
elseif hh > 0
    s = sprintf('%02dh %02dm %02ds', hh, mm, ss);
elseif mm > 0
    s = sprintf('%02dm %02ds', mm, ss);
else
    s = sprintf('%02ds', ss);
end
end
