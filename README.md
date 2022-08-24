# trusscalc_2D
MATLAB code that solves a plannar (2D) truss using the displacement (stiffness) method.

# Solution of a sample truss
A typical plannar truss is solved using the code. This truss is used as the default example in the code.


## Truss Geometry
<img src="https://github.com/FK-MAD/trusscalc_2D/blob/main/images/example%20truss.png" width="700">


## Boundary conditions

### Supports
Node 1 is pinned and node 9 is free to move only on the x-axis (x-roller).

### Loads
Nodes 3,5 and 7 are subjected to a force of 50 kN with direction in the negative y-axis.


## Results
      Node   x displacement [m]   y displacement [m]     x Force [N]     y Force [N]
         1                    0                    0      -7.713e-10         7.5e+04
         2              0.00025            -0.000303               0               0
         3             3.75e-05           -0.0005686               0          -5e+04
         4             0.000175           -0.0007237               0               0
         5             0.000125           -0.0007914               0          -5e+04
         6              7.5e-05           -0.0007237               0               0
         7            0.0002125           -0.0005686               0          -5e+04
         8            5.514e-19            -0.000303               0               0
         9              0.00025                    0               0         7.5e+04

-------------------------------------------------------------------------------------

    Member   Area [m^2]       E [Pa]       Force [N]          Strain     Stress [Pa]
         1     1.00e-02     2.00e+11      -1.061e+05      -5.303e-05      -1.061e+07
         2     1.00e-02     2.00e+11         7.5e+04        3.75e-05         7.5e+06
         3     1.00e-02     2.00e+11       1.061e+05       5.303e-05       1.061e+07
         4     1.00e-02     2.00e+11        -1.5e+05        -7.5e-05        -1.5e+07
         5     1.00e-02     2.00e+11      -3.536e+04      -1.768e-05      -3.536e+06
         6     1.00e-02     2.00e+11        1.75e+05        8.75e-05        1.75e+07
         7     1.00e-02     2.00e+11       3.536e+04       1.768e-05       3.536e+06
         8     1.00e-02     2.00e+11          -2e+05         -0.0001          -2e+07
         9     1.00e-02     2.00e+11       3.536e+04       1.768e-05       3.536e+06
        10     1.00e-02     2.00e+11        1.75e+05        8.75e-05        1.75e+07
        11     1.00e-02     2.00e+11      -3.536e+04      -1.768e-05      -3.536e+06
        12     1.00e-02     2.00e+11        -1.5e+05        -7.5e-05        -1.5e+07
        13     1.00e-02     2.00e+11       1.061e+05       5.303e-05       1.061e+07
        14     1.00e-02     2.00e+11         7.5e+04        3.75e-05         7.5e+06
        15     1.00e-02     2.00e+11      -1.061e+05      -5.303e-05      -1.061e+07

<img src="https://github.com/FK-MAD/trusscalc_2D/blob/main/images/example%20truss%20-%20displacement.png" width="700">

# License
This work is licensed under a
[Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International Public License][cc-by-nc-nd].

[![CC BY-NC-ND 4.0][cc-by-nc-nd-image]][cc-by-nc-nd]

[cc-by-nc-nd]: http://creativecommons.org/licenses/by-nc-nd/4.0/
[cc-by-nc-nd-image]: https://licensebuttons.net/l/by-nc-nd/4.0/88x31.png
[cc-by-nc-nd-shield]: https://img.shields.io/badge/License-CC%20BY--NC--ND%204.0-lightgrey.svg
