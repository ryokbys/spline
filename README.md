# Cubic Spline fitting

## Usage

1. prepare data file `dat.test` like following,
   ```
   4           ! num of data
   0.0  0.0    ! x1, y1
   1.0  2.0    ! x2, y2
   2.0  1.0    ! ...
   3.0  2.0    ! ...
   ```

2. Compile
   ```bash
   $ ifort -o spline cubic_spline.F lasubs.F
   ```

3. Run
   ```bash
   $ ./spline < dat.test
   ```
   Then you can get a list of fitted parameters.


## Note

* Results are written in standard output as following,
  ```
    0.0000000E+00  2.9333333E+00  0.0000000E+00 -9.3333333E-01
   -2.6000000E+00  1.0733333E+01 -7.8000000E+00  1.6666667E+00
    1.6600000E+01 -1.8066667E+01  6.6000000E+00 -7.3333333E-01
  ```
  (This output is the result of cubic spline fitting of `dat.test` above.)
  This output means,
  ```
    a0(1)  a1(1)  a2(1)  a3(1)
    a0(2)  a1(2)  a2(2)  a3(2)
    a0(3)  a1(3)  a2(3)  a3(3)
  ```
  The segment-1 which is in between data \\( (x1,y1) \\) and \\( (x2,y2) \\) is fitted to
  $$ y= a0(1) +a1(1)*x +a2(1)*x**2 +a3(1)*x**3. $$
  So the cubic polynomial form of segment-i is
  $$  y= sum_{j=0,3} aj(i)*x^j. $$

