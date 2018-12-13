function pde = Stokes_pipe
%% SINCOSDATA trigonometric  data for Poisson equation
%
% A simple model of colliding flow. The force f = 0, the velocity  u1 =
% 1-y^2, u2 = 0, p = -2x.
%
% Parabolic inflow boundary condition, natural outflow boundary condition.
% The constant in p is chosen such that the Neumann boundary condition
% du_1/dx - p = on x = 1.
%
% Reference: page 217 5.1.1 in Finite Elements and Fast Iterative Solvers with
% Applications in Incompressible Fluid Dynamics. by Howard C. Elman, David
% J. Silvester, and Andrew J. Wathen.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%'f',        @f : Load data (right hand side function)
%'exactu',   @exactu : Exact solution of u
%'g_D',      @g_D :  Dirichlet boundary condition
%'g_N',      @g_N :  Neumann boundary condition
%'Du',       @Du : Derivative of the exact solution u 
%'viscosity',@Diff1 : Diffusion coefficient tensor in case of scalar
%'Df',       @Df: Derivative of the right hand side function
%'pp',       @pp: Exact solution of p
%
%  Copyright (C)  Yujie LIU. Junping WANG. See COPYRIGHT.txt for details.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
    % load data (right hand side function)
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function rhs =  f(p)
    x = p(:,1); y = p(:,2);
    rhs(:,1) = x*0.0; %0.0 means no source term
    rhs(:,2) = x*0.0;
    rhs(:,3) = x*0.0;
    end
    % Dirichlet boundary condition
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function u =  g_D(p)
    u =  exactu(p);
    end
    % Neumann boundary condition
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function rhs =  g_N(p)
        x = p(:,1); y = p(:,2);
        rhs(:,1) = 2.*x./x;
        rhs(:,2) = x*0.0;
        rhs(:,3) = x*0.0;
    end
    % exact solution of u
    KnownSol= 1; % set this =0 if no exact solution is known
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function u = exactu(p)
    x = p(:,1); y = p(:,2);
    u(:,1) = 1-y.^2;
    u(:,2) = 0.*x;
    end
    
    % Derivative of the exact solution u
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function uprime =  Du(p)
    x = p(:,1); y = p(:,2);
    uprime(:,1,1) = 0.*x;
    uprime(:,1,2) = -2.*y;
    uprime(:,2,1) = 0.*x;
    uprime(:,2,2) = 0.*x;
    end
    % Diffusion coefficient tensor in case of scalar: This is the viscosity
    % for STOKES test
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function coef =  Diff1(p)
    x = p(:,1); y = p(:,2);
    coef = 0.*x+1;
    end
    % Derivative of the right hand side function
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function rhs =  Df(p)
    x = p(:,1); y = p(:,2);
    rhs(:,1,1) = 0.*x;
    rhs(:,1,2) = 0.*x;
    rhs(:,2,1) = 0.*x;
    rhs(:,2,2) = 0.*x;
    rhs(:,3,1) = 0.*x;
    rhs(:,3,2) = 0.*x;
    end
    % Exact solution of p
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function rhs =  pp(p)
    x = p(:,1); y = p(:,2);
    rhs =-2.*x;
    end
%
%pde = struct('f',@f,'exactu',@exactu,'g_D',@g_D,'Du',@Du, ...
%             'viscosity',@Diff1,'Df',@Df,'pp',@pp, 'KnownSol', KnownSol);
pde = struct('f',@f,'exactu',@exactu,'g_D',@g_D,'g_N',@g_N,'Du',@Du, ...
             'viscosity',@Diff1,'Df',@Df,'pp',@pp, 'KnownSol', KnownSol);
end