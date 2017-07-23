# Shadow Puppet Theatre Monte Carlo Ray-tracer

Creating a virtual shadow puppet theatre using monte-carlo ray-tracing for part of my EngD dissertation.

![puppet](https://user-images.githubusercontent.com/25514442/28500725-2df07684-6fc5-11e7-8c9d-2579553c2948.png)

![puppet](https://user-images.githubusercontent.com/25514442/28500759-dd4d913e-6fc5-11e7-8910-12b6df5b5549.png)

### Input:
 * quad for puppet
 * Bitmap as screen texture
 * Bitmap of scanned puppet as puppet texture
 
 ### Method
 This software can produce results for:
 * Rectangular light sources
 * Spherical light sources.
 
 ### Output:
 
 This software produces two models. The first for a small intense light source such a phone light and the other from a larger source such as a desk lamp.
 
 Also present is bounding boxes using a binary search tree, Monte-Carlo ray-tracing using adaptive sampling and multi-threading.
