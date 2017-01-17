function [convective] = setupConvective(E21, u, u_pres_orig, Ht02, h)


    % Set up txi vector which contains the vorticities of the
    % outer-oriented 0-cell
    xi = E21*u + u_pres_orig;
    txi = Ht02*xi;
    
    
    
    
    
    
    convective = [convection_u, convection_v];  