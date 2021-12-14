# GGP-for-Additive-Manufacturing

Based off from the recently published article in Ref[1], GGP is a MATLAB based code that has successfully integrated major geometric projection based topology optimizer, i.e., Geometric Projection (GP), MMC (Moving Morphable Components) and MNA (Moving Node Approach). This repository implements the Additive Layer Manufacturing (or 3D printing) constraints into GGP.

## How To cite us

@article{bhat2021some,
  title={On some applications of Generalized Geometric Projection to optimal 3D printing},
  author={Bhat, Krishnaraj Vilasraj and Capasso, Gabriele and Coniglio, Simone and Morlier, Joseph and Gogu, Christian},
  journal={Computers \& Graphics},
  year={2021},
  publisher={Elsevier}
}

## How To Use

Running the code is simple enough! Navigate to the directory in MATLAB file browser where the files have been downloaded, and use the following command in the MATLAB console to run structural topology optimizer:

```bash
GGP_AM(100,50,0.4,'MBB','MNA',8,3,10,-pi/2,-10)
```
1) The first two arguements dictate the resolution of the design space (100 x 50 elements in this case)
2) The third argument (0.4 in example) provides the maximum volume fraction constraint
3) The 4th argument mentions the boundary conditions ('MBB' in example)
4) The 5th argument specifies the projection method.
5) The 6th argument dictates the number of constraints
6) The 7th argument provides the maximum bridge length constraint
7) The last argument provides the initial longitudinal displacement of the design space.

Read below to know more about the arguments.

### Boundary Conditions

The 4th argument for the structural topology optimization is for the boundary condition, i.e., the problem set-up. Of course custom BCs can be used, but for reference, there are several example BCs already used namely:

1) 'MBB' : The classic Messerschmitt-Bolkow-Blohm beam.

2) 'L_Shape' : The L-Shaped Cantilever beam.

3) 'Short_Cantilever' : The simple short Cantilever beam.

### Method
The 5th argument for both structural and thermal boundary conditions is the method of solving. This can be changed to 3 methods:

1) 'GP' : Also known as Geometric Projection, first envisioned by Julian Norato (Ref[2]).

2) 'MMC' : Moving Morphable Components, by Xu Guo and team (Ref[3])

3) 'MNA' : Moving Node Approach, first conceptualised by J.T. Overvelde (Ref[4]), further developed in house at ISAE SUPAERO, Toulouse (France).

### Number of Constraints
This parameter in a whole integer that provides the user the freedom to apply the constraints from a chosen set of possibilities:
1) "1" : Applies only the volume fraction constraint.
2) "2" : Applies the volume fraction and overhang constraint (set at 45).
3) "3" : Applies the volume fraction, overhang and bridge length constraint. Note that the maximum bridge length constraint input is considered only when this option is used.

## Note:

1) While GP and MNA can be used, the MMC method is a work in progress.
2) If the user desires to read the output directly on the console, comment the lines 470 & 471.
3) The program automatically creates a folder with basic parameters of any run, and stores the density and component plot of selected iterations. The files and variables recorded also are stored here.
4) The user is adviced to run a smaller resolution if th intention is to mock execute the code, say 50 x 50. A full run with 100 x 100 elements takes well over 2 hours.

## References

[1. Coniglio, S., Morlier, J., Gogu, C., and Amargier, R., 2019.“Generalized geometry projection: A unified approach forgeometric feature based topology optimization”.](https://link.springer.com/article/10.1007/s11831-019-09362-8)

[2. Norato, J., Bell, B., and Tortorelli, D. A., 2015. “A geometry projection method for continuum-based topology optimization with discrete elements”.](https://www.sciencedirect.com/science/article/pii/S0045782515001711?casa_token=Xr892VegDc0AAAAA:89vzo5j0SLHYUh81j6ct9CI6nLxcAElsgHH-j3wqz5d1toX4X8BYiRwC3ZdPUg8Lu_Wyf3BtltM)

[3. Zhang, W., Yuan, J., Zhang, J., and Guo, X., 2016.“A new topology optimization approach based on moving morphable components (mmc) and the ersatz material model”.](https://idp.springer.com/authorize/casa?redirect_uri=https://link.springer.com/content/pdf/10.1007/s00158-015-1372-3.pdf&casa_token=iAKD3Y2P-30AAAAA:yMzRxgj07Jrk8lFPfZERQh7l05SX_PkJFCOmzNqBWRilfAOllY0mJ0dcDsOG7wX5qjq-66Ap8BkrqI2p)

[4. Overvelde, J. T., 2012. “The moving node approach intopology optimization”.](https://repository.tudelft.nl/islandora/object/uuid:86c056d8-f368-4239-893f-07ca3a22e112/datastream/OBJ1/download)
