%% Assignment for AE4134: CFD I
% Created by:
% Zhi-li Liu 4146557
% Jasper van Wensveen 4142179
%

function [u_pres] = setupupres(N, U_wall_top, U_wall_bot, U_wall_left, U_wall_right, ...
                                            V_wall_top, V_wall_bot, V_wall_left, V_wall_right, h)
% Start by creating arrays for each of the prescribed boundaries. For the u-velocities in the
% corners of the unit square it is assumed that velocities flowing tangent to a boundary belong
% to that specific boundary, and not to the one they are normal to.
u_pres_bot = ones(N+1,1)*U_wall_bot.*h';
u_pres_left = ones(N,1)*U_wall_left*h(1);
u_pres_right = ones(N,1)*U_wall_right*h(end);
u_pres_top = ones(N+1,1)*U_wall_top.*h';
v_pres_bot = ones(N,1)*V_wall_bot*h(1);
v_pres_left = ones(N+1,1)*V_wall_left.*h';
v_pres_right = ones(N+1,1)*V_wall_right.*h';
v_pres_top = ones(N,1)*V_wall_top*h(end);

% create u_press array alternating between left and right boundary
u_temp = [u_pres_left'; u_pres_right'];
u_pres_lr_mixed = u_temp(:);

% create v_press array alternating between left and right boundary
v_temp = [v_pres_left(2:end-1)'; v_pres_right(2:end-1)'];
v_pres_lr_mixed = v_temp(:);

% Construct the array with prescribed velocities                                      
u_pres = [ u_pres_bot;
           u_pres_lr_mixed;
           u_pres_top;
           v_pres_left(1);
           v_pres_bot;
           v_pres_right(1);
           v_pres_lr_mixed;
           v_pres_left(end);
           v_pres_top;
           v_pres_right(end) ];