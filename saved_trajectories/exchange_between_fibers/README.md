Help when visualising the trajectory:
+ Bonds are generated using cg_bond: ```cg_bonds -top diffusion_between_fibres.top```
+ The fibres are identified by ```index 0 to 6159``` and ```not index 0 to 6159```
+ The motor of interest is identified by ```index 10556 to 10583```
+ Another motor outside the fibre (that moved between the fibres during the more unstable equillibration phase) is identified by ```index 8288 to 8315```. This can simply be removed from the visualisation or colored according to the fibre it belongs to.
+ The most important movement (when the non-polar head moves from one fibre to the other) occurs between frames 2900 and 3500.