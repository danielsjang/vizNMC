%% Dawn-Dusk SSO RPO: Inspector at Delta-Inclination — ECI Animation
%
% Scenario:
%   Host:      circular dawn-dusk SSO (600 km, ~97 deg inclination)
%              RAAN = 90 deg so orbit normal points along +X_ECI (sun)
%              => sun perpendicular to orbital plane at equinox
%   Inspector: same RAAN, slightly different inclination (delta_inc_deg)
%              => out-of-plane NMC component; inspector swings toward/away
%                 from sun over each orbit
%              => near north pole: inspector tilts toward sun (exclusion)
%              => near south pole: inspector tilts away (clear or eclipse)
%
% With sun perpendicular to host orbit plane:
%   - Matched inclination inspector: LOS always 90 deg to sun => always clear
%   - Delta-inc inspector: periodic exclusion at one pole, clear at other
%
% ECI: X = vernal equinox (sun direction), Y = 90 deg east, Z = North pole

clear; clc; close all;

%% =========================================================================
%  USER PARAMETERS
% =========================================================================
alt_km             = 600;      % host circular altitude [km]
inclination_deg    = 97.0;     % host SSO inclination [deg]
delta_inc_deg      = 1.0;      % inspector inclination offset [deg]
                               % positive = inspector inclined more toward sun side
                               % try: 0.5, 1, 2, 5 deg for varying effect
RAAN_deg           = 90;       % dawn-dusk: orbit normal along +X_ECI (sun)
                               % do not change unless exploring other geometries
exclusion_half_deg = 60;       % optical exclusion cone half-angle [deg]
N_frames           = 240;      % frames per orbit
N_orbits           = 1;        % orbits to animate
pause_time         = 0.02;     % seconds between frames
ECI_exaggerate     = 12;       % relative motion scale for visibility
                               % (out-of-plane NMC is larger than in-plane
                               %  for small delta_inc, so less exaggeration needed)
save_gif           = true;    % set true to save gif (requires gif toolbox)
gif_filename       = 'DawnDusk_SSO_RPO.gif';
gif_fps            = 50;

%% =========================================================================
%  DERIVED PARAMETERS
% =========================================================================
R_earth      = 6371;
mu           = 398600.4418;
R_orbit      = R_earth + alt_km;
omega        = sqrt(mu / R_orbit^3);
T_orbit      = 2*pi / omega;

N_total = N_frames * N_orbits;
t_vec   = linspace(0, N_orbits * T_orbit, N_total);
t_min   = t_vec / 60;

rho_earth_deg      = asind(R_earth / R_orbit);
cos_eclipse_thresh = -sind(rho_earth_deg);
cos_excl           = cosd(exclusion_half_deg);

inc_host = deg2rad(inclination_deg);
inc_insp = deg2rad(inclination_deg + delta_inc_deg);   % inspector inclination
RAAN     = deg2rad(RAAN_deg);

% -------------------------------------------------------------------------
%  DAWN-DUSK RAAN VERIFICATION
%  N_hat = [sin(RAAN)*sin(i), -cos(RAAN)*sin(i), cos(i)]
%  For RAAN=90: N_hat = [sin(i), 0, cos(i)]
%  dot(N_hat, sun_ECI=[1,0,0]) = sin(i) ~ sin(97deg) ~ 0.993
%  => sun nearly along orbit normal. Not exactly 90 deg in-plane
%  because 97 deg inclination slightly offsets it, but close enough
%  for the dawn-dusk approximation.
% -------------------------------------------------------------------------
sun_ECI = [1; 0; 0];

N_hat_check = [sin(RAAN)*sin(inc_host); -cos(RAAN)*sin(inc_host); cos(inc_host)];
sun_dot_N   = dot(N_hat_check, sun_ECI);
sun_elev_from_plane = asind(sun_dot_N);   % elevation above orbital plane [deg]

fprintf('--- Dawn-Dusk SSO Parameters ---\n');
fprintf('  Altitude:            %d km\n',    alt_km);
fprintf('  Host inclination:    %.1f deg\n', inclination_deg);
fprintf('  Inspector delta-inc: %.2f deg\n', delta_inc_deg);
fprintf('  RAAN:                %.1f deg\n', RAAN_deg);
fprintf('  Period:              %.1f min\n', T_orbit/60);
fprintf('  Sun elev above plane: %.1f deg (90=perpendicular=ideal dawn-dusk)\n', ...
        sun_elev_from_plane);
fprintf('  Excl half-angle:     %.1f deg\n', exclusion_half_deg);

%% =========================================================================
%  PRECOMPUTE ORBITS
%
%  Host: circular, inc_host, RAAN
%  Inspector: circular, inc_insp, same RAAN, same SMA
%             => purely out-of-plane relative motion from delta_inc
%
%  CW relative motion for delta-inclination (no eccentricity difference):
%    R(t) = 0                                    [no radial offset]
%    T(t) = 0                                    [no along-track offset at t=0]
%    N(t) = R_orbit * delta_inc * sin(omega*t)   [out-of-plane oscillation]
%
%  This is the linearized result for a pure inclination difference with
%  nodes aligned. The inspector crosses the host's plane at t=0 (equator
%  crossing at ascending node) and reaches max N displacement at t=T/4.
%
%  NOTE: In practice a delta-inc also causes a small T drift unless
%  corrected, but for short arcs and small delta_inc this is negligible.
% =========================================================================
delta_inc    = deg2rad(delta_inc_deg);
A_N          = R_orbit * delta_inc;          % out-of-plane amplitude [km]
A_N_display  = A_N * ECI_exaggerate;        % exaggerated for ECI plot

fprintf('  Out-of-plane amplitude: %.2f km (true), %.0f km (displayed x%d)\n', ...
        A_N, A_N_display, ECI_exaggerate);

% Preallocate
host_ECI   = zeros(3, N_total);
insp_ECI   = zeros(3, N_total);
sun_R_unit = zeros(1, N_total);
sun_T_unit = zeros(1, N_total);
sun_N_unit = zeros(1, N_total);
nmc_R_rel  = zeros(1, N_total);
nmc_T_rel  = zeros(1, N_total);
nmc_N_rel  = zeros(1, N_total);

for k = 1:N_total
    u = omega * t_vec(k);   % host argument of latitude

    %% Host RTN basis in ECI
    R_hat = [ cos(RAAN)*cos(u) - sin(RAAN)*sin(u)*cos(inc_host);
              sin(RAAN)*cos(u) + cos(RAAN)*sin(u)*cos(inc_host);
              sin(u)*sin(inc_host) ];
    N_hat = [  sin(RAAN)*sin(inc_host);
              -cos(RAAN)*sin(inc_host);
               cos(inc_host) ];
    T_hat = cross(N_hat, R_hat);

    %% Sun in host RTN
    M = [R_hat'; T_hat'; N_hat'];
    s = M * sun_ECI;
    sun_R_unit(k) = s(1);
    sun_T_unit(k) = s(2);
    sun_N_unit(k) = s(3);

    %% Host ECI position
    host_ECI(:,k) = R_orbit * R_hat;

    %% Inspector relative position in host RTN
    %  Pure delta-inc: N oscillation only (linearized)
    nmc_R_rel(k) = 0;
    nmc_T_rel(k) = 0;
    nmc_N_rel(k) = A_N * sin(u);   % peaks at u=90 (north pole), trough at u=270 (south)

    %% Inspector absolute ECI (exaggerated N offset)
    rel_RTN       = [nmc_R_rel(k); nmc_T_rel(k); nmc_N_rel(k)];
    RTN2ECI       = [R_hat, T_hat, N_hat];
    insp_ECI(:,k) = host_ECI(:,k) + RTN2ECI * rel_RTN * ECI_exaggerate;
end

%% =========================================================================
%  ECLIPSE AND EXCLUSION FLAGS
%  Exclusion: angle between host->inspector LOS and sun < exclusion_half_deg
%  For dawn-dusk, sun is along N_hat. Inspector is offset in N.
%  => LOS is along N_hat direction => nearly aligned with sun near poles
%     => exclusion near north pole (inspector on sun side)
%     => clear or eclipse near south pole
% =========================================================================
in_eclipse   = false(1, N_total);
in_exclusion = false(1, N_total);

for k = 1:N_total
    in_eclipse(k) = (sun_R_unit(k) < cos_eclipse_thresh);

    if ~in_eclipse(k)
        pos_k = [nmc_R_rel(k); nmc_T_rel(k); nmc_N_rel(k)];
        nm    = norm(pos_k);
        if nm > 1e-6
            pn = pos_k / nm;
            sn = [sun_R_unit(k); sun_T_unit(k); sun_N_unit(k)];
            in_exclusion(k) = (dot(pn,sn) >= cos_excl);
        end
    end
end

ecl_pct  = 100*sum(in_eclipse)/N_total;
excl_pct = 100*sum(in_exclusion)/N_total;
fprintf('  Eclipse fraction:   %.1f%%\n', ecl_pct);
fprintf('  Exclusion fraction: %.1f%%\n', excl_pct);
fprintf('  Clear fraction:     %.1f%%\n', 100-ecl_pct-excl_pct);

%% =========================================================================
%  VISUAL EARTH RADIUS
%  Inspector is displaced in N (out of plane). Max displacement from Earth
%  center = sqrt(R_orbit^2 + A_N_display^2). Min along-track distance from
%  Earth center = R_orbit (when N=0). Scale Earth to clear both.
% =========================================================================
R_insp_min  = R_orbit - A_N_display * 0.5;   % conservative min clearance
R_earth_vis = min(R_earth, 0.72 * R_orbit);  % never bigger than true, stays clear

fprintf('  Visual Earth radius: %.0f km\n', R_earth_vis);

%% =========================================================================
%  FIGURE — ECI FRAME
% =========================================================================
fig = figure('Name','Dawn-Dusk SSO RPO — Delta-Inclination Inspector, ECI Frame', ...
             'Color',[0.06 0.06 0.10], 'Position',[100 80 600 450]);

ax = axes('Parent',fig, ...
          'Color',    [0.06 0.06 0.10], ...
          'XColor',   [0.65 0.65 0.65], ...
          'YColor',   [0.65 0.65 0.65], ...
          'ZColor',   [0.65 0.65 0.65], ...
          'GridColor',[0.22 0.22 0.22], ...
          'GridAlpha',0.5);
hold(ax,'on'); grid(ax,'on'); axis(ax,'equal');

lim = R_orbit * 1.55;
xlim(ax,[-lim lim]); ylim(ax,[-lim lim]); zlim(ax,[-lim lim]);
xlabel(ax,'X_{ECI}  [km]','Color',[0.80 0.80 0.80],'FontSize',11);
ylabel(ax,'Y_{ECI}  [km]','Color',[0.80 0.80 0.80],'FontSize',11);
zlabel(ax,'Z_{ECI}  [km]','Color',[0.80 0.80 0.80],'FontSize',11);
% title(ax, {'Dawn-Dusk SSO RPO — \DeltaInclination Inspector, ECI Frame', ...
%            sprintf('%d km, %.0f° Incl, RAAN=%.0f°  |  \\Deltai = %.2f°  |  Offset ×%d  |  Equinox', ...
%            alt_km, inclination_deg, RAAN_deg, delta_inc_deg, ECI_exaggerate), ...
%            'Blue = Clear  |  Red = Solar Exclusion  |  Dark Blue = Eclipse'}, ...
%       'Color',[1.0 1.0 0.65],'FontSize',11);
view(ax, 30, 22);   % lower elevation to better show out-of-plane motion

%% -------------------------------------------------------------------------
%  STATIC ELEMENTS
%% -------------------------------------------------------------------------

% Earth sphere
[xe,ye,ze] = sphere(60);
surf(ax, R_earth_vis*xe, R_earth_vis*ye, R_earth_vis*ze, ...
    'FaceColor',[0.12 0.28 0.55],'EdgeColor','none','FaceAlpha',0.92);

% Equatorial ring
th_eq = linspace(0,2*pi,200);
plot3(ax, R_earth_vis*cos(th_eq), R_earth_vis*sin(th_eq), zeros(1,200), ...
    'Color',[0.40 0.55 0.80 0.5],'LineWidth',1.0);

% North/South pole markers
plot3(ax,0,0, R_earth_vis,'w.','MarkerSize',8);
plot3(ax,0,0,-R_earth_vis,'w.','MarkerSize',8);
text(ax,0,0, R_earth_vis+lim*0.05,'N','Color',[0.7 0.7 0.7],'FontSize',9, ...
    'HorizontalAlignment','center');
text(ax,0,0,-R_earth_vis-lim*0.07,'S','Color',[0.7 0.7 0.7],'FontSize',9, ...
    'HorizontalAlignment','center');

% Full host orbit trace (faint)
plot3(ax, host_ECI(1,:), host_ECI(2,:), host_ECI(3,:), '--', ...
    'Color',[0.45 0.65 0.90 0.18],'LineWidth',0.8);

% Full inspector orbit trace (faint)
plot3(ax, insp_ECI(1,:), insp_ECI(2,:), insp_ECI(3,:), '--', ...
    'Color',[0.85 0.55 0.20 0.15],'LineWidth',0.8);

% Sun arrow — along +X_ECI, perpendicular to dawn-dusk orbit plane
sun_len = lim * 0.88;
quiver3(ax,0,0,0, sun_len,0,0, 0, ...
    'Color',[1.0 0.85 0.0],'LineWidth',2.8,'MaxHeadSize',0.25);
text(ax, sun_len*1.03, 0, lim*0.04,'SUN', ...
    'Color',[1.0 0.85 0.0],'FontSize',12,'FontWeight','bold');

% Dawn/dusk labels — sun perpendicular to orbit, so ±Y is dawn/dusk
text(ax, 0,  lim*0.88, 0,'DUSK','Color',[0.65 0.50 0.80],'FontSize',9, ...
    'FontWeight','bold','HorizontalAlignment','center');
text(ax, 0, -lim*0.88, 0,'DAWN','Color',[0.65 0.50 0.80],'FontSize',9, ...
    'FontWeight','bold','HorizontalAlignment','center');

% Note: for dawn-dusk orbit, ascending node is at +Y (RAAN=90)
% Sun is along +X (into the page from orbit perspective)
% orbit goes over north pole (+Z) and south pole (-Z)

% ECI axis triad
al = R_earth_vis * 0.50;
o  = [-lim*0.82, -lim*0.82, -lim*0.82];
quiver3(ax,o(1),o(2),o(3), al,0,0, 0,'r','LineWidth',1.5,'MaxHeadSize',1.2);
quiver3(ax,o(1),o(2),o(3), 0,al,0, 0,'g','LineWidth',1.5,'MaxHeadSize',1.2);
quiver3(ax,o(1),o(2),o(3), 0,0,al, 0,'b','LineWidth',1.5,'MaxHeadSize',1.2);
text(ax,o(1)+al*1.3,o(2),      o(3),      'X','Color','red',  'FontSize',9,'FontWeight','bold');
text(ax,o(1),      o(2)+al*1.3,o(3),      'Y','Color','green','FontSize',9,'FontWeight','bold');
text(ax,o(1),      o(2),      o(3)+al*1.3,'Z','Color','cyan', 'FontSize',9,'FontWeight','bold');

%% -------------------------------------------------------------------------
%  DYNAMIC HANDLES
%% -------------------------------------------------------------------------
h_host   = plot3(ax, host_ECI(1,1), host_ECI(2,1), host_ECI(3,1), 'w^', ...
    'MarkerSize',11,'MarkerFaceColor','white','MarkerEdgeColor','white');
h_insp   = plot3(ax, insp_ECI(1,1), insp_ECI(2,1), insp_ECI(3,1), 'o', ...
    'MarkerSize',9,'MarkerFaceColor',[0.3 0.7 1.0],'MarkerEdgeColor','white');
h_los    = plot3(ax, [host_ECI(1,1) insp_ECI(1,1)], ...
                     [host_ECI(2,1) insp_ECI(2,1)], ...
                     [host_ECI(3,1) insp_ECI(3,1)], ...
    'Color',[0.3 0.7 1.0],'LineWidth',1.6);

N_trail  = 80;
h_htrail = plot3(ax, host_ECI(1,1), host_ECI(2,1), host_ECI(3,1), ...
    'Color',[0.60 0.78 1.00 0.45],'LineWidth',1.3);
h_itrail = plot3(ax, insp_ECI(1,1), insp_ECI(2,1), insp_ECI(3,1), ...
    'Color',[1.00 0.62 0.25 0.45],'LineWidth',1.3);

h_status = text(ax,-lim*0.95, lim*0.72, lim*0.92,'', ...
    'FontSize',13,'FontWeight','bold');
h_time   = text(ax,-lim*0.95, lim*0.72, lim*0.80,'', ...
    'Color',[0.85 0.85 0.85],'FontSize',10);
h_orbit  = text(ax,-lim*0.95, lim*0.72, lim*0.70,'', ...
    'Color',[0.80 0.80 0.50],'FontSize',10);
h_Ndisp  = text(ax,-lim*0.95, lim*0.72, lim*0.60,'', ...
    'Color',[0.70 0.70 0.70],'FontSize',9);

%% =========================================================================
%  ANIMATION LOOP
% =========================================================================
for k = 1:N_total

    eclipsed = in_eclipse(k);
    excluded = in_exclusion(k);

    if eclipsed
        c_insp     = [0.15 0.20 0.55];
        c_los      = [0.20 0.25 0.60];
        status_str = '🌑  ECLIPSE — No Optical';
        status_col = [0.40 0.45 0.90];
    elseif excluded
        c_insp     = [1.00 0.25 0.20];
        c_los      = [1.00 0.25 0.20];
        status_str = '⚠  SOLAR EXCLUSION';
        status_col = [1.00 0.35 0.30];
    else
        c_insp     = [0.30 0.70 1.00];
        c_los      = [0.30 0.70 1.00];
        status_str = '✓  LOS Clear';
        status_col = [0.35 1.00 0.50];
    end

    set(h_host,'XData',host_ECI(1,k),'YData',host_ECI(2,k),'ZData',host_ECI(3,k));
    set(h_insp,'XData',insp_ECI(1,k),'YData',insp_ECI(2,k),'ZData',insp_ECI(3,k), ...
               'MarkerFaceColor',c_insp);
    set(h_los, 'XData',[host_ECI(1,k) insp_ECI(1,k)], ...
               'YData',[host_ECI(2,k) insp_ECI(2,k)], ...
               'ZData',[host_ECI(3,k) insp_ECI(3,k)],'Color',c_los);

    k0 = max(1,k-N_trail);
    set(h_htrail,'XData',host_ECI(1,k0:k),'YData',host_ECI(2,k0:k),'ZData',host_ECI(3,k0:k));
    set(h_itrail,'XData',insp_ECI(1,k0:k),'YData',insp_ECI(2,k0:k),'ZData',insp_ECI(3,k0:k));

    set(h_status,'String',status_str,'Color',status_col);

    elapsed_min = t_min(k);
    orb_num     = floor(elapsed_min/(T_orbit/60))+1;
    orb_phase   = mod(elapsed_min, T_orbit/60);
    set(h_time, 'String',sprintf('T = %.1f min  (phase %.1f / %.1f min)', ...
                 elapsed_min, orb_phase, T_orbit/60));
    set(h_orbit,'String',sprintf('Orbit %d / %d', orb_num, N_orbits));
    set(h_Ndisp,'String',sprintf('N offset: %.1f km (true)  |  Δi = %.2f°', ...
                 nmc_N_rel(k), delta_inc_deg));

    drawnow;
    pause(pause_time);

    % GIF export (requires gif toolbox)
    if save_gif
        if k == 1
            pause();
            gif(gif_filename,'DelayTime',1/gif_fps);
        else
            gif;
        end
    end
end

hold(ax,'off');
