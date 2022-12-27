# Image Compression

The idea behind this project is to implement a barely simplified version of the image codification algorithm for `.jpg` archives, which is based on an essential property of the Fourier Transform: transforming signals from real life usually leads to low frequencies with little contribution from high frequencies. This notion introduces the following idea: given a signal we can consider it's Fourier Transform and discard higher frequencies. This allows us to store only a small number of frequencies, the compressed signal. To retrieve the signal you can complete the transformed vector with zeros and do the inverse transform.

This algorithm uses the Discrete Cosine Transform, mainly due to the fact that it converts Real numbers into Real numbers.

## Preparation

This simplification of the algorithm asumes the image dimensions are divisibles by 16. Therefore, the folowing function fills the image's margins with black to ensure that both dimensions are divisible by 16.

## Using the right format

The algorithm uses the YCbCr decomposition. Therefore, first we must convert the image to the correct format. We also separate the three channels to operate with them properly. Moreover, due to the human eye being muvh more sensitive to the luminosity than to the intensity of a color, we can reduce the matrices Cb and Cr by creating smaller matrices in which every pixel is the average of 4 neighbouring pixels in the original matrix. Finally, we'll center all of the matrices coefficients on 0. The values will range between -128 and 128.

The inverse process is also implemented, since it will be used later to retrieve the compressed image.

## Using the Discrete Cosine Transform

The following step consists in thinking of each matrix as a group of `8×8` submatrices and applying the discrete fourier transform to each.

Similarly to the previous step, the inverse process is also implemented, since it will be used later to retrieve the compressed image.
