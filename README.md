# Ciwinska_Nature_2024


`mammary.f` is a Fortran code that models the clonal dynamics of the mouse mammary gland epithelium.

The dynamics are described in more detail in Ciwinska et al. main text and Supplementary Information. In short, the mammary gland ductal epithelium is modeled as a 1D chain with each cell labeled with a unique barcode at day zero. To model turnover during rounds of the estrous cycle, a domain of cells is chosen at random and removed. Cells on either side of the domain then duplicate and expand to replenish the domain. To avoid commensurability effects, we further introduce a degree of cell loss and replacement. Details of the parameters can be found in the Fortran code.

## Running the Code

The code can be run on a Unix platform with the following command line:

```sh
gfortran mammary.f
./a.out


The outputs are provided as clsize.dat and cldist.dat, with the former containing the average statistics as a function of estrous cycle number and the latter containing the clone size distributions. 
