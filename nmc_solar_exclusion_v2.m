%% NMC Solar Exclusion Visualization v2 — Physically Correct, GEO RTN Frame
%
% Host is at RTN origin (GEO circular orbit).
% Inspector follows CW (Clohessy-Wiltshire) natural motion circumnavigation.
% Sun vector rotates in the R-T plane at the same angular rate as the GEO orbit (once/sidereal day).
%
% PHASE CONVENTION:
%   phase_offset_deg = 0   -> Unfavorable: NMC east (+T) when sun is east (+T)
%   phase_offset_deg = 180 -> Favorable:   NMC east (+T) when sun is west (-T)
%
% RTN:  R = radial (zenith), T = tangential (east/velocity dir), N = orbit normal

clear; clc; close all;

%% -------------------------------------------------------------------------
%  USER PARAMETERS
% -------------------------------------------------------------------------
phase_offset_deg   = 0;       % 0 = unfavorable, 180 = favorable [deg]
exclusion_half_deg = 60;      % optical sensor solar exclusion half-angle [deg]
semi_minor_R       = 5;       % NMC radial semi-minor axis [km]
N_frames           = 240;     % total animation frames (one full GEO day)
pause_time         = 0.02;    % seconds between frames

%% -------------------------------------------------------------------------
%  DERIVED PARAMETERS
% -------------------------------------------------------------------------
semi_major_T = 2 * semi_minor_R;   % CW gives 2:1 T:R naturally [km]
T_GEO        = 86164;              % GEO sidereal period [sec]
omega        = 2*pi / T_GEO;       % mean motion [rad/s]
t_vec        = linspace(0, T_GEO, N_frames);  % time vector [sec]
t_hrs        = t_vec / 3600;       % time in hours for display

phase_offset_rad = deg2rad(phase_offset_deg);
cos_excl         = cosd(exclusion_half_deg);

%% -------------------------------------------------------------------------
%  CW NATURAL MOTION CIRCUMNAVIGATION (no drift, pure NMC)
%
%  Initial conditions for a CW NMC starting at theta_0:
%    x0  =  a_R * sin(theta_0)          [radial, R]
%    y0  =  a_T * cos(theta_0)          [tangential, T]  (note: 2*a_R = a_T)
%    xd0 =  a_R * omega * cos(theta_0)
%    yd0 = -a_T * omega * sin(theta_0)  (= -2*a_R*omega*sin for 2:1)
%
%  CW solution for NMC (zero drift, C=0):
%    R(t) =  a_R * sin(omega*t + theta_0)
%    T(t) =  2*a_R * cos(omega*t + theta_0)   [= a_T * cos(...)]
%    N(t) =  0
%
%  theta_0 = 0 -> starts at T+ (east) at t=0, consistent with unfavorable phase
% -------------------------------------------------------------------------
theta_0_nmc = 0;   % NMC starts east (+T) at t=0
nmc_R_pos   =  semi_minor_R * sin(omega * t_vec + theta_0_nmc);
nmc_T_pos   =  semi_major_T * cos(omega * t_vec + theta_0_nmc);
nmc_N_pos   =  zeros(1, N_frames);

%% -------------------------------------------------------------------------
%  SUN VECTOR — rotates in R-T plane at GEO rate
%  Sun also starts at +T (east) at t=0 for unfavorable phase (phase_offset=0)
%  sun_angle = omega*t + phase_offset
%    T-component = cos(sun_angle)  -> +T at t=0 when phase_offset=0
%    R-component = sin(sun_angle)  -> small radial excursion as it rotates
% -------------------------------------------------------------------------
sun_angle  = omega * t_vec + phase_offset_rad;
sun_T      = cos(sun_angle);
sun_R      = sin(sun_angle);
sun_N      = zeros(1, N_frames);
sun_arrow_len = 13;   % visual length [km]

%% -------------------------------------------------------------------------
%  FIGURE SETUP
% -------------------------------------------------------------------------
fig = figure('Name', 'NMC Solar Exclusion v2 — RTN Frame', ...
             'Color', [0.07 0.07 0.10], ...
             'Position', [80 80 600 400]);

ax = axes('Parent', fig, ...
          'Color',      [0.07 0.07 0.10], ...
          'XColor',     [0.65 0.65 0.65], ...
          'YColor',     [0.65 0.65 0.65], ...
          'ZColor',     [0.65 0.65 0.65], ...
          'GridColor',  [0.25 0.25 0.25], ...
          'GridAlpha',  0.5);
hold(ax, 'on');
grid(ax, 'on');
axis(ax, 'equal');

lim = 15;
limz = lim/2;
xlim(ax, [-lim lim]);
ylim(ax, [-lim lim]);
zlim(ax, [-limz limz]);
xlabel(ax, 'Radial  R  [km]',     'Color', [0.80 0.80 0.80], 'FontSize', 11);
ylabel(ax, 'Tangential  T  [km]', 'Color', [0.80 0.80 0.80], 'FontSize', 11);
zlabel(ax, 'Normal  N  [km]',     'Color', [0.80 0.80 0.80], 'FontSize', 11);

if phase_offset_deg == 0
    phase_str = 'Unfavorable Phase (NMC east = Sun east)';
else
    phase_str = sprintf('Phase offset = %d°', phase_offset_deg);
end
% title(ax, {'GEO RPO — CW Natural Motion Circumnavigation, RTN Frame', phase_str}, ...
      % 'Color', [1.0 1.0 0.65], 'FontSize', 13);
view(ax, -30, 22);

%% -------------------------------------------------------------------------
%  STATIC ELEMENTS
% -------------------------------------------------------------------------

% Full NMC orbit trail (faint dashed)
plot3(ax, nmc_R_pos, nmc_T_pos, nmc_N_pos, '--', ...
    'Color', [0.35 0.55 0.85 0.35], 'LineWidth', 1.0);

% Host at origin
plot3(ax, 0, 0, 0, 'w^', 'MarkerSize', 13, 'MarkerFaceColor', 'white');
text(ax, 0.6, 0.6, 1.4, 'HOST', ...
    'Color', 'white', 'FontSize', 9, 'FontWeight', 'bold');

% Earth arrow (always -R in RTN)
earth_len = 11;
quiver3(ax, 0,0,0, -earth_len,0,0, 0, ...
    'Color', [0.25 0.75 0.35], 'LineWidth', 2.2, 'MaxHeadSize', 0.5);
text(ax, -earth_len*1.15, 0, 2, 'EARTH', ...
    'Color', [0.25 0.75 0.35], 'FontSize', 10, 'FontWeight', 'bold');

% Small RTN reference triad (bottom-left corner)
o = [-lim+1.5, -lim+1.5, -limz+1.5];
alen = 3.5;
quiver3(ax, o(1),o(2),o(3), alen,0,0, 0, 'r', 'LineWidth',1.5,'MaxHeadSize',1.0);
quiver3(ax, o(1),o(2),o(3), 0,alen,0, 0, 'g', 'LineWidth',1.5,'MaxHeadSize',1.0);
quiver3(ax, o(1),o(2),o(3), 0,0,alen, 0, 'b', 'LineWidth',1.5,'MaxHeadSize',1.0);
text(ax, o(1)+alen+0.5, o(2),       o(3),       'R', 'Color','red',   'FontSize',9,'FontWeight','bold');
text(ax, o(1),          o(2)+alen+0.5, o(3),     'T', 'Color','green', 'FontSize',9,'FontWeight','bold');
text(ax, o(1),          o(2),       o(3)+alen+0.5,'N', 'Color','cyan',  'FontSize',9,'FontWeight','bold');

%% -------------------------------------------------------------------------
%  DYNAMIC HANDLES (updated each frame)
% -------------------------------------------------------------------------

% Sun arrow (dynamic — rotates)
sv = sun_arrow_len * [sun_R(1); sun_T(1); sun_N(1)];
h_sun = quiver3(ax, 0,0,0, sv(1),sv(2),sv(3), 0, ...
    'Color', [1.0 0.85 0.0], 'LineWidth', 2.8, 'MaxHeadSize', 0.45);
h_sun_lbl = text(ax, sv(1)*1.15, sv(2)*1.15, sv(3)+1.5, 'SUN ', ...
    'Color', [1.0 0.85 0.0], 'FontSize', 11, 'FontWeight', 'bold');

% Inspector (NMC object)
h_nmc = plot3(ax, nmc_R_pos(1), nmc_T_pos(1), 0, 'o', ...
    'MarkerSize', 11, ...
    'MarkerFaceColor', [0.3 0.7 1.0], ...
    'MarkerEdgeColor', 'white', ...
    'LineWidth', 1.2);

% LOS line host -> NMC
h_los = plot3(ax, [0 nmc_R_pos(1)], [0 nmc_T_pos(1)], [0 0], ...
    'Color', [0.3 0.7 1.0], 'LineWidth', 1.8);

% Trailing history of NMC (last N_trail points)
N_trail = 40;
h_trail = plot3(ax, nmc_R_pos(1), nmc_T_pos(1), 0, ...
    'Color', [0.3 0.7 1.0 0.4], 'LineWidth', 1.2);

% Status text (top-left)
h_status = text(ax, -lim+0.4,  lim-1.0,  limz-0.3, '', ...
    'FontSize', 12, 'FontWeight', 'bold');

% Time text
h_time = text(ax, -lim+0.4,  lim-1.0,  limz-2.2, '', ...
    'Color', [0.85 0.85 0.85], 'FontSize', 11);

% Sun-NMC angle readout
h_angle = text(ax, -lim+0.4,  lim-1.0,  limz-4.0, '', ...
    'Color', [0.75 0.75 0.75], 'FontSize', 10);

%% -------------------------------------------------------------------------
%  ANIMATION LOOP
% -------------------------------------------------------------------------
for k = 1 : N_frames

    % --- NMC position ---
    pos = [nmc_R_pos(k); nmc_T_pos(k); nmc_N_pos(k)];

    % --- Sun direction unit vector ---
    sun_unit = [sun_R(k); sun_T(k); sun_N(k)];   % already unit length

    % --- LOS unit vector from host to NMC ---
    pos_norm = pos / (norm(pos) + 1e-12);

    % --- Solar exclusion check ---
    cos_angle    = dot(pos_norm, sun_unit);
    angle_deg    = acosd(max(-1, min(1, cos_angle)));
    in_exclusion = (cos_angle >= cos_excl);

    % Colors
    if in_exclusion
        nmc_color = [1.00 0.25 0.20];
        los_color = [1.00 0.25 0.20];
        set(h_status, 'String', ' ⚠  SOLAR EXCLUSION  ⚠ ', 'Color', [1.0 0.30 0.25]);
    else
        nmc_color = [0.30 0.70 1.00];
        los_color = [0.30 0.70 1.00];
        set(h_status, 'String', '✓  LOS Clear', 'Color', [0.35 1.0 0.50]);
    end

    % --- Update NMC marker ---
    set(h_nmc, 'XData', nmc_R_pos(k), ...
               'YData', nmc_T_pos(k), ...
               'ZData', nmc_N_pos(k), ...
               'MarkerFaceColor', nmc_color);

    % --- Update LOS line ---
    set(h_los, 'XData', [0, nmc_R_pos(k)], ...
               'YData', [0, nmc_T_pos(k)], ...
               'ZData', [0, nmc_N_pos(k)], ...
               'Color', los_color);

    % --- Update trail ---
    k0 = max(1, k - N_trail);
    set(h_trail, 'XData', nmc_R_pos(k0:k), ...
                 'YData', nmc_T_pos(k0:k), ...
                 'ZData', nmc_N_pos(k0:k));

    % --- Update sun arrow ---
    sv = sun_arrow_len * sun_unit;
    set(h_sun, 'UData', sv(1), 'VData', sv(2), 'WData', sv(3));
    set(h_sun_lbl, 'Position', [sv(1)*1.10, sv(2)*1.10, sv(3)+1.0]);

    % --- Update text overlays ---
    hrs = floor(t_hrs(k));
    mins = floor((t_hrs(k) - hrs) * 60);
    set(h_time,  'String', sprintf('Time:  %02d h %02d m  /  24 h', hrs, mins));
    set(h_angle, 'String', sprintf('LOS–Sun angle:  %.1f°  (excl < %d°)', ...
                                    angle_deg, exclusion_half_deg));

    drawnow;
    pause(pause_time);
    if k == 1
        pause
        gif('GEO_badNMC.gif','DelayTime',1/50)
    else
        gif
    end
end

hold(ax, 'off');
