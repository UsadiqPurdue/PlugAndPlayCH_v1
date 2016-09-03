function [sino,ANG_D,ANG_R]=SetSinoGeom(max,min,offset,rr,step)

[ANG_D,ANG_R,step]=GetUniformAngles(max,min,offset,step);

% Initialize Sinogram's geometry


% Tilt Settings
n_t=length(ANG_D);

sino.n_theta=n_t;
sino.theta_0=max.*pi/180;
sino.delta_theta=-step;sino.delta_t=1;

% delta_t and delta_theta remain the same = 1°

sino.cosine=cos(ANG_R);
sino.sine=sin(ANG_R);

% Displacement Settings

sino.n_t=rr;
sino.t_0=-1*rr/2; % Center of the detector

end