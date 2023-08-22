# 2D frequency dependant resistance and inductance calculation method

This method is based on analytical solution under quasimagneto static field for round condcutors.

The possible impacts from magnetic core and air gaps is considered.

This code is built in MATLAB 2019b

**Example**

There are 2 conductors with 1mm radius. They locate at (0,0) and (2.2e-3,0), unit m, respectively.

To calculate the impedance at 1MHz with opposite 1A in conductors, and the order is set to 5, commend is listed as follow

x = [0;2.2e-3]; y = [0;0]; a = 1e-3*ones(2,1); f = 1e6; I = [1,-1];

Z = MultiConMatrix_L(x,y,a,'freq',f,'I',I,'Nord',5);

Result Z is 

0.08793 + 0.65456i

0.08793 + 0.65456i

Result from COMSOL is

0.08836+0.65416i

0.08836+0.65416i

**Author:**

Tianming Luo, TU Delft, HV group

https://www.tudelft.nl/ewi/over-de-faculteit/afdelingen/electrical-sustainable-energy/high-voltage-technologies/research/design-method-for-high-power-medium-frequency-transformers