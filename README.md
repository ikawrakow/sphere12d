# Sampling random points in a 12D sphere

## 1. Why?

A while ago I was giving lectures at a course about Monte Carlo (MC) simulations in Medical Physics using [EGSnrc](https://nrc-cnrc.github.io/EGSnrc/). One of my lectures was on Variance Reduction Techniques in the context of radiation transport simulations, and I was looking for a good (but simple) example to demonstrate to the students how one can increase the efficiency of MC simulations by orders of magnitude by some modest amount of cleverness. One classic introductory example for MC is the computation of $\pi$. One knows that the area of a circle of radius $R$ is $A_c = \pi R^2$. The area of a square of half-size $R$ is $A_s = 4 R^2$, so the ratio of the two is $A_c/A_s = \pi/4$. Hence, if we draw $N$ random points within a square with half-size $R$ (simply `x = 2*rndm() - 1, y = 2*rndm() - 1`, where `rndm()` is a function returning random values uniformly distributed in `0...1`), and count the number of points $N_{\rm in}$ inside the circle, we have $\pi = 4 p$, where $p = N_{\rm in}/N$ is the estimated probability for the point to be inside the circle. We can even get the statistical uncertainty of the estimate, i.e., $\pi = 4 p \pm \sqrt{p (1 - p)/(N-1)}$. OK, that's nice and simple, but I cannot impress the students with that. So I came up with the idea to look at drawing random points in a high-dimensional sphere, say 12D. Why? Let's see.

## 2. Simple (rejection method)

In analogy to the circle, the simplest implementation to pick a random point in a 12D unit sphere is
1. Pick 12 random numbers $x_i$ uniformly distributed in `[-1, 1]`.
2. Compute $r^2 = \sum x_i^2$.
3. If $r^2 \le 1$ deliver $\{x_i\}$
4. Goto step 1

Let's see how this goes. The code is in `simple12d.cpp`. I'm using a [Xorshift](https://en.wikipedia.org/wiki/Xorshift) random number generator (RNG), specifically `xoshiro256+`, which is currently my favorite RNG (passes all known tests, very fast). The code is in `randomGenerator.h` and `randomGenerator.cpp`. So
```
cmake --build .
./bin/simple12d
It took 306201633 attempts to sample 100000 points in a 12d sphere
Sampling efficiency: 0.000326582
<r^2> = 0.856875
Run time: 5675.12 ms
```
This is on a Ryzen-7950X CPU. We need to sample ~3,000 random 12d points in a 12d cube to find one that is inside the 12d sphere! This is in stark contrast to 2d (circle, fraction of points inside is $\pi/4$) or a sphere (3d, fraction of points inside is $\pi/6$).   Did you know that the volume of a $K$-dimensional unit sphere is such a minuscule fraction of the volume of the $K$-dimensional cube in which it is inscribed when $K$ is large? We can consult [Wikipedia](https://en.wikipedia.org/wiki/Volume_of_an_n-ball): volume of a 12d unit sphere is $V_{s,12} = \pi^6/720$, volume of a 12d cube with unit half-size is $V_{c,12} = 2^{12}$, so $V_{s,12}/V_{c,12} = (\pi/4)^6/720 \approx 0.000326$, just like our MC simulation predicts. With that result, I can now go and impress the students by devising a method that is 500+ times faster

## 3. Smart

We need to do some math. The probability distribution function (pdf) for points in a 12d unit sphere is

$${\rm d}F(x_1, x_2, \cdots, x_{12}) = \Theta \left(1 - \sum_{i=1}^{12} x_i^2\right) \prod_{i=1}^{12} {\rm d}x_i,~~~~ -1 \le x_i \le 1$$

where $\Theta$ is the Heaviside step function. Let's go to polar coordinates for pairs of points, i.e.
```math
\begin{eqnarray}
x_{2 i + 0} & = r_i \cos \phi_i \\
x_{2 i + 1} & = r_i \sin \phi_i
\end{eqnarray}
```
which results in

$${\rm d} F(r_1, \phi_1, \cdots r_6, \phi_6) = \Theta \left(1 - \sum_{i=1}^{6} r_i^2\right) \prod_{i=1}^{6} r_i {\rm d}r_i {\rm d}\phi_i,~~~~ 0 \le r_i \le 1,~ 0 \le \phi_i \le 2 \pi$$

Now apply the trivial variable transformation $y_i = r_i^2$ and then change variables again according to
```math
\begin{eqnarray}
z_1 & = & y_1\\
z_2 & = & y_1 + y_2\\
z_3 & = & y_1 + y_2 + y_3\\
& \cdots &\\
z_6 & = & y_1 + y_2 + y_3 + y_4 + y_5 + y_6
\end{eqnarray}
```
which results in

$${\rm d}F(z_1, \phi_1, \cdots, z_6, \phi_6) = (1/2)^6~\Theta \left(1 - z_6\right) \prod_{i=1}^{6} {\rm d}z_i {\rm d}\phi_i,~~~~ z_1 \le z_2 \le z_3 \le z_4 \le z_5 \le z6$$

I.e., all we need to do to get the $\{z_i}$ is to pick 6 random numbers in `[0, 1]`, sort them in ascending order, and assign the first to $z_1$, the second to $z_2$, etc. To check our math we can compute the volume of the sphere by integrating over all variables. The integration over the azimuthal angles $\phi_i$ is trivial and gives a factor of $2 \pi^6$, so we have

$$V_{s, 12} = \pi^6 \int_0^1 {\rm d}z_6 \int_{0}^{z_6} {\rm d}z_5 \int_{0}^{z_5} {\rm d}z_4 \cdots \int_{0}^{z_2} {\rm d}z_1 = \pi^6/6! = \pi^6/720$$

as expected from the Wikipedia article. The sampling algorithm is then as follows
1. Pick 6 random numbers in `[0, 1]`
2. Sort in ascending order and assign the lowest to $z_1$, second lowest to $z_2$, etc.
3. Set $r_i = \sqrt{z_i - z_{i-1}}$ where we use the convention $z_0 = 0$.
4. Pick 6 azimuthal angles randomly in $\[0, 2 \pi\]$
5. Deliver $x_{2i + 0} = r_i \cos \phi_i,~x_{2i + 1} = r_i \sin \phi_i,~i = 1 \cdots 6$

without the need for any rejections. This is implemented in `smart12d.cpp` and we get
```
./bin/smart12d
It took 100000 attempts to sample 100000 points in a 12d sphere
Sampling efficiency: 1
<r^2> = 0.85702
Run time: 10.627 ms
```
Nice. `5675 / 10.6 = 537` times faster than the simple method.

I'm computing and reporting the average $r^2$ if the points inside the sphere to prevent the compiler from optimizing out the calls to the sampling function in the time measurement loop. But is `<r^2> = 0.85702` actually correct? As $r^2 = \sum r_i^2 = \sum y_i \equiv z_6$, we can compute $<r^2>$ analytically from

$$< r^2 > = \pi^6/V_{s, 12} = \int_0^1 z_6 {\rm d}z_6 \int_{0}^{z_6} {\rm d}z_5 \int_{0}^{z_5} {\rm d}z_4 \cdots \int_{0}^{z_2} {\rm d}z_1 = 6/7 \approx 0.8571$$

so yes, the result is correct (within statistical uncertainty).

Obviously one can easily generalize the code here to sample random points in any $K$ dimensional sphere as long as $K$ is even. When $K$ is odd it becomes slightly more complicated because one has a spare coordinate where one cannot use polar coordinates. I leave the extension to odd $K$ for another day.

## 4. Circles and 3D spheres

We started by talking about picking random points inside a circle and using a rejection method. Obviously one can go to polar coordinates and have a rejection free method:
1. Pick random $r$ in `[0...1]` using $r {\rm d}r$ as pdf. I.e., $r = \sqrt{\eta}$ where $\eta$ is a random number uniformly distributed in `[0...1]`.
2. Pick a random angle in $\[0, 2 \pi\]$
3. Deliver point $(r \cos \phi, r \sin \phi)$

 But is this faster? `circle.cpp` has several methods for picking random points in a (unit) circle. We get
 ```
./bin/circle

===================== Rejection(1000000 points)
<r2> = 0.500039
Time: 8.54295 ms

===================== Direct1(1000000 points)
<r2> = 0.499589
Time: 27.7229 ms

===================== Direct2(1000000 points)
<r2> = 0.50024
Time: 9.45063 ms

===================== Direct3(1000000 points)
<r2> = 0.499968
Time: 10.5975 ms
```

Clearly no. Rejection is still fastest. `Direct1` is the above algorithm, so more than 3 times slower (despite hardware implementations for `sqrt` and `cos/sin`). `Direct2` uses a different method for obtaining $\cos \phi$ and $\sin \phi$ directly, without the need to evaluate trigonometric functions (see the `randomAzimuth()` function in `circle.cpp`). This is much faster than `Direct1`, but still slower than rejection. `Direct3` uses a trick to avoid the evaluation of `sqrt`: one picks 2 random numbers in `[0, 1]` and uses the larger of the two for $r$ in step 1 of the above algorithm (left as exercise to prove that this gives the correct pdf for $r$). In the old days `max(rndm(), rndm()` used to be faster than `sqrt(rndm())`, but the `sqrt` hardware implementation on the Ryzen is aparently fast enough to beat `max(rndm(), rndm())`.

What about 3D?
```
./bin/sphere

===================== Rejection(1000000 points)
<r2> = 0.600087
Time: 23.8875 ms

===================== Direct1(1000000 points)
<r2> = 0.600181
Time: 20.7947 ms

===================== Direct2(1000000 points)
<r2> = 0.600106
Time: 13.8537 ms
```
Here rejection becomes slower. The fastest method (o my CPU) is this:
1. Set `r = max(rndm(), max(rndm(), rndm()))` (maximum of 3 random numbers uniformly distributed in `[0, 1]`. The resulting pdf for `r` is $r^2 {\rm d}r$
2. Set `cost = 2*rndm()-1, sint = sqrt(1 - cost*cost)` (cosine and sine of the polar angle)
3. Pick `cphi, sphi` using `randomAzimuth()`
4. Deliver 3D point `(r*sint*cphi, r*sint*sphi, r*cost)`
