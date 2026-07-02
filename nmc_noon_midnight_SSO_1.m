%% Worst-Case Optical RPO: Noon-Midnight SSO — ECI Animation
%
% Scenario:
%   Host:      circular noon-midnight SSO (600 km, ~97 deg inclination)
%              RAAN chosen so sun (+X_ECI, equinox) lies exactly in the
%              orbital plane => sun along ±T in RTN (worst case for optics)
%   Inspector: same inclination/RAAN, small eccentricity
%              perigee toward sun (+T) => closest to host at noon (excluded)
%              apogee away from sun (-T) => farthest at midnight (eclipse)
%
% Result: close approach = solar exclusion, far = eclipse.
%         Only brief clear windows at dawn/dusk quadrants.
%
% ECI: X = vernal equinox, Y = 90 deg east equatorial, Z = North pole
%
% NOTE: Earth sphere is scaled DOWN visually so exaggerated inspector
%       orbit traces always clear the displayed sphere.

clear; clc; close all;

%% =========================================================================
%  USER PARAMETERS
% =========================================================================
alt_km             = 600;      % host circular altitude [km]
inclination_deg    = 97.0;     % SSO inclination [deg]
ecc_inspector      = 0.004;    % inspector eccentricity (small, visible in ECI)
omega_peri_deg     = 0;       % arg of perigee: 90 deg => perigee at
                               % top of orbit = noon/sun side
exclusion_half_deg = 60;       % optical exclusion cone half-angle [deg]
N_frames           = 240;      % frames per orbit
N_orbits           = 1;        % orbits to animate
pause_time         = 0.02;     % seconds between frames
ECI_exaggerate     = 60;       % relative motion scale factor for visibility

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

inc        = deg2rad(inclination_deg);
omega_peri = deg2rad(omega_peri_deg);

% -------------------------------------------------------------------------
%  COMPUTE RAAN SO SUN LIES IN THE ORBITAL PLANE
%
%  Orbit normal in ECI:  N_hat = [sin(RAAN)*sin(i), -cos(RAAN)*sin(i), cos(i)]
%  Sun in ECI:           sun_ECI = [1, 0, 0]
%
%  Sun in orbital plane  <=>  N_hat . sun_ECI = 0
%    => sin(RAAN)*sin(i) = 0
%    => RAAN = 0 or 180 deg  (for i ~= 0)
%
%  RAAN = 0:   ascending node at +X_ECI (vernal equinox direction, same as sun)
%              => sun is at the ascending node, along +T at start
%  RAAN = 180: ascending node at -X_ECI
%              => sun is at the descending node
%
%  For noon-midnight: we want sun along +T (noon) when host is at the
%  ascending node crossing. RAAN=0 achieves this cleanly.
% -------------------------------------------------------------------------
RAAN_deg = 0;
RAAN     = deg2rad(RAAN_deg);

fprintf('--- Noon-Midnight SSO ---\n');
fprintf('  Altitude:         %d km\n',    alt_km);
fprintf('  Inclination:      %.1f deg\n', inclination_deg);
fprintf('  RAAN (computed):  %.1f deg\n', RAAN_deg);
fprintf('  Inspector ecc:    %.4f\n',     ecc_inspector);
fprintf('  Arg of perigee:   %.1f deg\n', omega_peri_deg);
fprintf('  Period:           %.1f min\n', T_orbit/60);
fprintf('  NMC radial amp:   %.1f km\n',  R_orbit*ecc_inspector);
fprintf('  ECI exaggeration: x%d\n',      ECI_exaggerate);

%% =========================================================================
%  SUN VECTOR — fixed in ECI at equinox
% =========================================================================
sun_ECI = [1; 0; 0];

%% =========================================================================
%  PRECOMPUTE ORBITS AND SUN IN RTN
% =========================================================================
a_R_nmc  = R_orbit * ecc_inspector;
a_T_nmc  = 2 * a_R_nmc;
phi_nmc  = omega_peri + pi/2;   % phase: perigee at omega_peri

host_ECI  = zeros(3, N_total);
insp_ECI  = zeros(3, N_total);
sun_R_unit = zeros(1, N_total);
sun_T_unit = zeros(1, N_total);
sun_N_unit = zeros(1, N_total);
nmc_R_rel  = zeros(1, N_total);
nmc_T_rel  = zeros(1, N_total);

for k = 1:N_total
    u = omega * t_vec(k);   % argument of latitude

    % ECI basis vectors for RTN
    R_hat = [ cos(RAAN)*cos(u) - sin(RAAN)*sin(u)*cos(inc);
              sin(RAAN)*cos(u) + cos(RAAN)*sin(u)*cos(inc);
              sin(u)*sin(inc) ];
    N_hat = [  sin(RAAN)*sin(inc);
              -cos(RAAN)*sin(inc);
               cos(inc) ];
    T_hat = cross(N_hat, R_hat);

    % Sun in RTN
    M = [R_hat'; T_hat'; N_hat'];
    s = M * sun_ECI;
    sun_R_unit(k) = s(1);
    sun_T_unit(k) = s(2);
    sun_N_unit(k) = s(3);

    % Host ECI
    host_ECI(:,k) = R_orbit * R_hat;

    % NMC relative position (RTN)
    nmc_R_rel(k) =  a_R_nmc * sin(u + phi_nmc);
    nmc_T_rel(k) =  a_T_nmc * cos(u + phi_nmc);

    % Inspector ECI (exaggerated offset)
    rel_RTN      = [nmc_R_rel(k); nmc_T_rel(k); 0];
    RTN2ECI      = [R_hat, T_hat, N_hat];
    insp_ECI(:,k) = host_ECI(:,k) + RTN2ECI * rel_RTN * ECI_exaggerate;
end

%% =========================================================================
%  ECLIPSE AND EXCLUSION FLAGS
% =========================================================================
in_eclipse   = false(1, N_total);
in_exclusion = false(1, N_total);

for k = 1:N_total
    in_eclipse(k) = (sun_R_unit(k) < cos_eclipse_thresh);
    if ~in_eclipse(k)
        pos_k = [nmc_R_rel(k); nmc_T_rel(k); 0];
        pn    = pos_k / (norm(pos_k) + 1e-12);
        sn    = [sun_R_unit(k); sun_T_unit(k); sun_N_unit(k)];
        in_exclusion(k) = (dot(pn,sn) >= cos_excl);
    end
end

ecl_pct  = 100*sum(in_eclipse)/N_total;
excl_pct = 100*sum(in_exclusion)/N_total;
fprintf('  Eclipse fraction:   %.1f%%\n', ecl_pct);
fprintf('  Exclusion fraction: %.1f%%\n', excl_pct);
fprintf('  Clear fraction:     %.1f%%\n', 100-ecl_pct-excl_pct);

%% =========================================================================
%  VISUAL EARTH RADIUS
%  Scale Earth sphere so it is comfortably inside the inspector orbit trace.
%  Inspector ECI positions are at R_orbit +/- ECI_exaggerate*a_R_nmc from
%  center. Minimum inspector radius ~ R_orbit - ECI_exaggerate*a_R_nmc.
%  Set visual Earth radius to 80% of that minimum.
% =========================================================================
R_insp_min  = R_orbit - ECI_exaggerate * a_R_nmc;
R_earth_vis = 0.75 * R_insp_min;

fprintf('  True Earth radius:   %.0f km\n', R_earth);
fprintf('  Visual Earth radius: %.0f km\n', R_earth_vis);
fprintf('  Min inspector dist:  %.0f km\n', R_insp_min);

%% =========================================================================
%  FIGURE — ECI FRAME
% =========================================================================
fig = figure('Name','ECI Frame — Noon-Midnight SSO Worst Case', ...
             'Color',[0.06 0.06 0.10], 'Position',[100 80 600 450]);

ax = axes('Parent',fig, ...
          'Color',    [0.06 0.06 0.10], ...
          'XColor',   [0.65 0.65 0.65], ...
          'YColor',   [0.65 0.65 0.65], ...
          'ZColor',   [0.65 0.65 0.65], ...
          'GridColor',[0.22 0.22 0.22], ...
          'GridAlpha',0.5);
hold(ax,'on'); grid(ax,'on'); axis(ax,'equal');

lim = R_orbit * 1.5;
xlim(ax,[-lim lim]); ylim(ax,[-lim lim]); zlim(ax,[-lim lim]);
xlabel(ax,'X_{ECI}  [km]','Color',[0.80 0.80 0.80],'FontSize',11);
ylabel(ax,'Y_{ECI}  [km]','Color',[0.80 0.80 0.80],'FontSize',11);
zlabel(ax,'Z_{ECI}  [km]','Color',[0.80 0.80 0.80],'FontSize',11);
% title(ax, {'Worst-Case Optical RPO — Noon-Midnight SSO, ECI Frame', ...
%            sprintf('%d km Circular, %.0f° Incl  |  Inspector offset ×%d  |  Equinox', ...
%            alt_km, inclination_deg, ECI_exaggerate), ...
%            'Blue = Clear  |  Red = Solar Exclusion  |  Dark Blue = Eclipse'}, ...
%       'Color',[1.0 1.0 0.65],'FontSize',11);
view(ax, 25, 30);

%% -------------------------------------------------------------------------
%  STATIC ELEMENTS
%% -------------------------------------------------------------------------

% Earth sphere — visually scaled down
[xe,ye,ze] = sphere(60);
surf(ax, R_earth_vis*xe, R_earth_vis*ye, R_earth_vis*ze, ...
    'FaceColor',[0.12 0.28 0.55], ...
    'EdgeColor','none', ...
    'FaceAlpha', 0.92);
% Subtle equatorial ring on sphere
th_eq = linspace(0,2*pi,200);
plot3(ax, R_earth_vis*cos(th_eq), R_earth_vis*sin(th_eq), zeros(1,200), ...
    'Color',[0.40 0.55 0.80 0.6],'LineWidth',1.0);
% North pole dot
plot3(ax,0,0,R_earth_vis,'w.','MarkerSize',8);
text(ax,0,0,R_earth_vis+lim*0.04,'N','Color',[0.7 0.7 0.7],'FontSize',8, ...
    'HorizontalAlignment','center');

% Full host orbit trace (faint)
plot3(ax, host_ECI(1,:), host_ECI(2,:), host_ECI(3,:), '--', ...
    'Color',[0.45 0.65 0.90 0.20],'LineWidth',0.8);

% Full inspector orbit trace (faint)
plot3(ax, insp_ECI(1,:), insp_ECI(2,:), insp_ECI(3,:), '--', ...
    'Color',[0.85 0.55 0.20 0.15],'LineWidth',0.8);

% Sun arrow — static, in +X_ECI, in orbital plane (RAAN=0, so orbital
% plane contains X axis — visually confirm sun lies in orbit plane)
sun_len = lim * 0.88;
quiver3(ax,0,0,0, sun_len,0,0, 0, ...
    'Color',[1.0 0.85 0.0],'LineWidth',2.8,'MaxHeadSize',0.25);
text(ax, sun_len*1.04, 0, lim*0.04, 'SUN', ...
    'Color',[1.0 0.85 0.0],'FontSize',12,'FontWeight','bold');

% "Midnight" label on anti-sun side
text(ax,-sun_len*0.88, 0, lim*0.04,'MIDNIGHT', ...
    'Color',[0.45 0.45 0.75],'FontSize',9,'FontWeight','bold', ...
    'HorizontalAlignment','center');

% ECI axis triad
al  = R_earth_vis * 0.55;
o   = [-lim*0.82, -lim*0.82, -lim*0.82];
quiver3(ax,o(1),o(2),o(3), al,0,0, 0,'r','LineWidth',1.5,'MaxHeadSize',1.2);
quiver3(ax,o(1),o(2),o(3), 0,al,0, 0,'g','LineWidth',1.5,'MaxHeadSize',1.2);
quiver3(ax,o(1),o(2),o(3), 0,0,al, 0,'b','LineWidth',1.5,'MaxHeadSize',1.2);
text(ax,o(1)+al*1.2,o(2),      o(3),      'X','Color','red',  'FontSize',9,'FontWeight','bold');
text(ax,o(1),      o(2)+al*1.2,o(3),      'Y','Color','green','FontSize',9,'FontWeight','bold');
text(ax,o(1),      o(2),      o(3)+al*1.2,'Z','Color','cyan', 'FontSize',9,'FontWeight','bold');

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

% Text overlays
h_status = text(ax,-lim*0.95, lim*0.72, lim*0.92,'', ...
    'FontSize',13,'FontWeight','bold');
h_time   = text(ax,-lim*0.95, lim*0.72, lim*0.80,'', ...
    'Color',[0.85 0.85 0.85],'FontSize',10);
h_orbit  = text(ax,-lim*0.95, lim*0.72, lim*0.70,'', ...
    'Color',[0.80 0.80 0.50],'FontSize',10);

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

    % Host marker
    set(h_host,'XData',host_ECI(1,k),'YData',host_ECI(2,k),'ZData',host_ECI(3,k));

    % Inspector marker + LOS
    set(h_insp,'XData',insp_ECI(1,k),'YData',insp_ECI(2,k),'ZData',insp_ECI(3,k), ...
               'MarkerFaceColor',c_insp);
    set(h_los, 'XData',[host_ECI(1,k) insp_ECI(1,k)], ...
               'YData',[host_ECI(2,k) insp_ECI(2,k)], ...
               'ZData',[host_ECI(3,k) insp_ECI(3,k)], ...
               'Color',c_los);

    % Trails
    k0 = max(1,k-N_trail);
    set(h_htrail,'XData',host_ECI(1,k0:k), ...
                 'YData',host_ECI(2,k0:k), ...
                 'ZData',host_ECI(3,k0:k));
    set(h_itrail,'XData',insp_ECI(1,k0:k), ...
                 'YData',insp_ECI(2,k0:k), ...
                 'ZData',insp_ECI(3,k0:k));

    % Text
    set(h_status,'String',status_str,'Color',status_col);

    elapsed_min = t_min(k);
    orb_num     = floor(elapsed_min/(T_orbit/60))+1;
    orb_phase   = mod(elapsed_min, T_orbit/60);
    set(h_time, 'String', sprintf('T = %.1f min  (phase %.1f / %.1f min)', ...
                 elapsed_min, orb_phase, T_orbit/60));
    set(h_orbit,'String', sprintf('Orbit %d / %d', orb_num, N_orbits));

    drawnow;
    pause(pause_time);
    if k == 1
        pause
        gif('LEOSSO_badNMC.gif','DelayTime',1/50)
    else
        gif
    end
end

hold(ax,'off');
