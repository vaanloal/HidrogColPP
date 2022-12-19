function [position, isterminal, direction] = quenchingEvent(t, x, Tfiltro)

    position(1) = Tfiltro - x(9); % The value that we want to be zero
    isterminal = 1; % Halt integration
    direction = 0; % The zero can be approached from either direction
end
