function result =  dgintcubature(field, cub)
	% Note: cub.W => cubature weights with Jacobian scalings.
	% cub.V => interpolation matrix to the cubature points.
    result = sum(sum(cub.W.*(cub.V*field)));
end
