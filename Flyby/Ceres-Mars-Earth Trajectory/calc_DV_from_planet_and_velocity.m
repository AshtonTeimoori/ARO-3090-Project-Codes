function delta_V = calc_DV_from_planet_and_velocity(vvec_planet, mu, v_inf, rp)

% Calculate circular orbit velocity
v_circ = sqrt(mu/rp);

% Calculate energy for hyperbolic trajectory
energy_hyp = (norm(v_inf)^2)/2;

% Calculate the velocity at perigee for hyperbolic trajectory
v_hyp = sqrt(2*(mu/rp+energy_hyp));

% Calculate the delta V
delta_V = v_hyp-v_circ;

end

