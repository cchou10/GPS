
%function to rotate coordinates to local vertical east north
function [V_diff E_diff N_diff]=ven(latOBS,longOBS,ECEFdiffx, ECEFdiffy, ECEFdiffZ)
V_diff = cos(latOBS)*(ECEFdiffx*cos(longOBS) + ...
		     ECEFdiffy*sin(longOBS)) + ECEFdiffZ*sin(latOBS);
	E_diff = ECEFdiffy*cos(longOBS) - ECEFdiffx*sin(longOBS);
	N_diff = ECEFdiffZ*cos(latOBS) - sin(latOBS)*...
		     (ECEFdiffx*cos(longOBS) + ECEFdiffy*sin(longOBS));

end