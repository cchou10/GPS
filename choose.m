% choose.m
%
%  Copyright (c) 2008 Paul M. Kinter, Jr., and 
%  Copyright (c) 2012 Mark L. Psiaki.  All rights reserved. 
%
% calculate 'a choose b': number of combinations of size 'b'
% for a group of 'a' numbers
%
% input: 'a' and 'b'
%				[ a, b]
%
function combos = choose(a,b)
	numer = factorial(a);
	denom = factorial(a-b)*factorial(b);
	combos = numer/denom;
return;