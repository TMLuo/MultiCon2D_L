# 2D frequency dependant resistance and inductance calculation method

This method is based on analytical solution under quasimagneto static field for round condcutors.

The possible impacts from magnetic core and air gaps is considered.

This code is built in MATLAB 2019b

**Example**

There are 2 conductors with 1mm radius. They locate at (0,0) and (2e-3,0), unit m, respectively.

To calculate the impedance at 1MHz with opposite 1A in conductors, commend is listed as follow

x = [0;2e-3]; y = [0;0]; a = 1e-3*ones(2,1); f = 1e6; I = [1,-1];

Z = MultiConMatrix_L(x,y,a,'freq',f,'I',I);

Result Z is 
0.0879314016907529 + 0.654564457663242i
0.0879314016907557 + 0.654564457663227i

Result from COMSOL
0.08835774737291231+0.6541584654416708i
0.08835774737291231+0.6541584654416708i

**Author:**

Tianming Luo, TU Delft, HV group