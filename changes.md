#	Wed. 3/2-21
edited minor definitions in headerfiles, essentially added "class" to some 
input variables to ensure the code still worked with them, as well as to 
have a more uniform description.

Currently assuming any non-variable factor is 1 -> h-bar, mass, etc.
only including alpha and r.

Restricting currently to 1 dimension and 1 particle. Can expand later.

change in Wavefunc/simplegaussian 
	compute derivative -> analytical solution for 1 particle in 1 dimension
	added. commented out exponent because they are expensive to calculate
	and I do not believe it will be necessary wrt. using the result in local 
	energy. It simply removes the necessity for factoring it away in the 
	local energy calculation.

Implemented local energy in "harmonicOscillator" double derivative and
evaluation of wavefunction in "simpleGaussian".
