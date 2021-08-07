# FORM--One-loop-integral-decompositon-
A FORM program to decompose four- leg loop integrals into Passarino-Veltman tensor integrals (A,B,C,D,... functions)
1. We use FORM-  a symbolic manipulation system as the tool to transform a four-leg loop integral into the sum of tensor integrals. So that we can use another numerical tool, e.g. COLLIER library, to solve the integrals in a particular problem. 
2.  As the very first version of the program, we limit the problem for the loop which just contains one kind of particle, i.e. the masses on four propagation lines are the same.
(Indeed, we have first applied the program to solve the Light- by- light scattering problem at leading order in QED,  where the loop is a closed fermion loop).
The output is represented in the standard form to be processed by COLLIER later.
3. For more information about FORM (manuals, tutorials,...). Visit website: [https://www.nikhef.nl/~form/](url). That is especially suitable for solving problems in particle physics.
4. COLLIER library is a Fortran- based program to give numerical results to the expressions which contain the Passarino-Veltman tensor integrals in a standard form.
This library can be applied up to 7-point functions with multiple options for the users. For more details about it, visit the website [https://collier.hepforge.org/](url)
