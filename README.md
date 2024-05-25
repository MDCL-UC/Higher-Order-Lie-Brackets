# Higher-Order-Lie-Brackets
This repo corresponds to project titled as "Higher-order Lie bracket approximation and averaging of control-affine systems with application to extremum seeking." Specifically, this repo consists of codes for Example 3 and Example 4 corresponding to the mentioned paper. Following is the description of each file in the folder.

For Example 3 (Or Example 4), params.m consists of essential parameters. Find_nu.m provides the nu that are needed for the example. Specifically, the function saves nu_val2.mat and nu_val3.mat that are loaded in the example.

Find_Beta.m provides beta needed for the example. It consists of functions beta2_fun() and beta3_fun(). The control inputs (u1,u2,u3,u4) must be checked within those functions. It saves beta_val.mat. VectorFieldsExample3.mlx( or VectorFieldsExample4.mlx) provides the symbolic expression for the Lie brackets.

Example3.m (or Example 4.m) consists of the main code that loads the mat file, simulates the system and runs the code
